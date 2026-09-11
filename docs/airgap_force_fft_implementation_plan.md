# ギャップ電磁力 2D FFT 実装計画

## 1. 目的

モータ電磁界解析から出力されるギャップ磁束密度 `Br`, `Btheta` を読み込み、
Maxwell 応力から以下を一貫して評価できる Python 処理系を実装する。

- 径方向電磁応力 `Fr`（実装上は `sigma_r` [Pa]）
- 接線方向電磁応力 `Ftheta`（実装上は `sigma_t` [Pa]）
- トルク
- 径方向合力 `Fx`, `Fy`
- UMP (Unbalanced Magnetic Pull)
- 空間次数 × 時間周波数の 2D FFT
- 空間次数 × ロータ回転次数の 2D FFT
- 支配的な電磁加振力モードの抽出
- 将来の振動解析・応力解析への入力データ生成

初期実装は NumPy / pandas / Matplotlib を使用する。
FFT コアは NumPy に限定し、依存関係を小さくする。

---

## 2. 背景と設計上の考え方

トルクリプルと径方向電磁加振力は別の指標である。

ギャップ磁束密度を

- `Br(theta, t)` : 径方向磁束密度
- `Btheta(theta, t)` : 周方向磁束密度

とすると、空気中の Maxwell 応力は

```text
sigma_r     = (Br^2 - Btheta^2) / (2 mu0)
sigma_theta = Br * Btheta / mu0
```

で評価する。

ここで、

- `sigma_theta` の空間 0 次成分は全周積分するとトルクになる。
- `sigma_r` の空間 1 次成分は正味の径方向力、すなわち UMP に対応する。
- `sigma_r` の空間 2 次成分は主として楕円変形を励起する。
- `sigma_r` の空間 3 次以上はより高次の円環変形モードを励起する。

したがって `Fr` と `Ftheta` は同じ 2D FFT エンジンで処理する。

---

## 3. 初期実装の入力仕様

### 3.1 推奨 CSV 形式

long-form CSV を標準とする。

```csv
time,theta_deg,Br,Btheta
0.000000,0.0,...,...
0.000000,1.0,...,...
0.000000,2.0,...,...
...
0.000100,0.0,...,...
```

必須情報:

- ステップ軸
  - `time` [s] または
  - `rotor_angle_deg` [deg]
- `theta_deg`
  - ステータ固定座標系でのギャップ周方向位置 [mechanical deg]
- `Br` [T]
- `Btheta` [T]

別途コマンドラインまたは設定ファイルから、

- ギャップ評価半径 `R` [m]
- 積厚 `L` [m]

を与える。

### 3.2 初期版の制約

初期版では以下を要求する。

1. `theta` は 0～360 deg の全周データ。
2. 360 deg 点は 0 deg と重複させない。
3. `theta` は等角度刻み。
4. `time` を使う場合は等時間刻み。
5. `rotor_angle_deg` を使う場合は等機械角刻み。
6. すべてのステップで同じ `theta` 点数を持つ。
7. NaN / 欠損点を許さない。

不等間隔データの再サンプリングは Phase 2 以降とする。

---

## 4. 座標系と符号規約

実装前に以下を明示的に固定する。

- `theta` 正方向
- ロータ回転正方向
- `Br` 正方向
- `Btheta` 正方向
- `Ftheta` 正方向
- 正トルク方向

この規約は README とテストデータに必ず残す。

2D FFT では符号付きの

- 空間次数 `m`
- 時間周波数 `f`

を保持する。

実数場では共役なピークが現れるが、初期版では正負を潰さない。
これにより進行波・後退波の区別を後から行える。

---

## 5. Maxwell 応力の計算

```text
mu0 = 4 pi 1e-7

sigma_r =
    (Br^2 - Btheta^2) / (2 mu0)

sigma_t =
    Br * Btheta / mu0
```

単位は Pa。

変数名として `Fr`, `Ftheta` を使う場合でも、
内部では力 [N] と応力 [Pa] を混同しないため

- `sigma_r`
- `sigma_t`

を推奨する。

---

## 6. 全周積分量

ギャップ円筒面の微小面積を

```text
dA = R * L * dtheta
```

とする。

### 6.1 トルク

```text
T(t) = R^2 * L * integral sigma_t(theta,t) dtheta
```

### 6.2 径方向合力

```text
Fx(t) =
    R * L * integral sigma_r(theta,t) cos(theta) dtheta

Fy(t) =
    R * L * integral sigma_r(theta,t) sin(theta) dtheta
```

### 6.3 UMP

```text
UMP(t) = sqrt(Fx^2 + Fy^2)
```

これらは FFT とは独立に直接積分して出力する。
FFT 実装の検証にも使用する。

---

## 7. 2D FFT

配列形状を

```text
field[step, theta]
```

と統一する。

対象 field:

- `sigma_r`
- `sigma_t`
- 将来的には `Br`
- 将来的には `Btheta`

NumPy `fft2` を使用し、

```text
C = fft2(field) / (Nstep * Ntheta)
```

とする。

`fftshift` を用いて符号付き軸を中央配置する。

### 7.1 空間軸

全周 360 deg データの場合、

```text
m = ..., -3, -2, -1, 0, 1, 2, 3, ...
```

を空間次数とする。

### 7.2 時間軸

`step=time` の場合:

```text
f [Hz]
```

### 7.3 ロータ角軸

`step=rotor_angle_mech` の場合:

```text
n [cycles / mechanical revolution]
```

すなわち rotor order とする。

ロボット用途のような変速運転では、
後述の角度領域 order tracking を優先する。

---

## 8. FFT 振幅の定義

初期版では `fft2 / (Nstep*Ntheta)` による
「複素 Fourier 係数」を正とする。

これは実装上最も曖昧さが少ない。

実数 cosine 波の物理振幅を表示する場合は、
共役ペアを考慮した `2*abs(C)` が必要になる場合がある。

したがって出力では以下を分ける。

- `complex_coefficient`
- `coefficient_abs`
- 将来追加: `cosine_amplitude`

DC、Nyquist、共役ペアの扱いを明示しないまま
一律に 2 倍しない。

---

## 9. 出力仕様

### 9.1 積分量

`integrated_metrics.csv`

```text
step
Fx_N
Fy_N
UMP_N
Torque_Nm
```

### 9.2 2D スペクトル

- `Fr_spectrum.npz`
- `Ftheta_spectrum.npz`

保存内容:

```text
step_frequency_or_order
spatial_order
complex_coefficient
```

複素データは CSV より NPZ を推奨する。

### 9.3 支配モード

- `Fr_dominant_modes.csv`
- `Ftheta_dominant_modes.csv`

列:

```text
step_frequency_or_order
spatial_order_m
complex_real
complex_imag
coefficient_abs
```

### 9.4 グラフ

初期版:

- `Fr_2d_fft.png`
- `Ftheta_2d_fft.png`

横軸:
- frequency [Hz] または rotor order

縦軸:
- spatial order `m`

色:
- `20 log10(|C|)`

---

## 10. 重要な検証項目

### Test A: 単一進行波

人工データ

```text
x(theta,t) =
    A cos(m0 theta - 2 pi f0 t)
```

を生成する。

期待結果:

- `(m0, -f0)`
- `(-m0, +f0)`

に共役ピーク。

符号は NumPy DFT 規約に合わせてテストで固定する。

### Test B: 一様接線応力

```text
sigma_t(theta,t) = const
```

期待結果:

- 空間 0 次のみ。
- 全周積分トルクと `m=0` Fourier 成分が一致。

### Test C: 空間 1 次径方向応力

```text
sigma_r(theta) = A cos(theta)
```

期待結果:

- `m=1` のみ。
- `Fx != 0`
- `Fy ~= 0`
- FFT の `m=1` と直接積分結果が整合。

### Test D: 空間 2 次径方向応力

```text
sigma_r(theta) = A cos(2 theta)
```

期待結果:

- `m=2` が存在。
- `Fx ~= 0`
- `Fy ~= 0`

すなわち楕円変形力は存在しても UMP は発生しないことを確認。

### Test E: FEM トルクとの比較

電磁界ソルバーが出力するトルクと

```text
T_from_gap =
    R^2 L integral sigma_t dtheta
```

を比較する。

差が大きい場合は

- ギャップ評価半径
- Br/Btheta の補間位置
- Maxwell stress の評価面
- 符号
- 単位
- 2D/3D の積厚扱い

を確認する。

---

## 11. 推奨モジュール構成

プロトタイプ確認後、Codex で以下に分割する。

```text
airgap_fft/
    __init__.py
    io.py
    stress.py
    metrics.py
    spectrum.py
    plotting.py
    config.py
    cli.py

tests/
    test_io.py
    test_stress.py
    test_metrics.py
    test_spectrum.py
    test_traveling_wave.py

examples/
    synthetic/
    gl80/
```

責務:

### `io.py`

- CSV reader
- 列名マッピング
- グリッド検査
- 単位変換

### `stress.py`

- Maxwell stress

### `metrics.py`

- torque
- Fx
- Fy
- UMP
- spatial mode amplitude

### `spectrum.py`

- 1D spatial FFT
- 1D temporal FFT
- 2D FFT
- mode extraction
- signed order handling

### `plotting.py`

- heatmap
- mode spectrum
- torque/UMP histories

### `cli.py`

- コマンドライン実行

---

## 12. 実装フェーズ

### Phase 0: 実データ形式確認

実際のモータ解析出力について確認する。

- 1 ファイルに全 step が入るか
- step ごとに別ファイルか
- `theta` は固定子座標かロータ座標か
- `Br`, `Btheta` の定義
- 0/360 deg 重複
- 全周か周期セクタか
- 時刻かロータ位置か
- 積厚、評価半径をどこから取得するか

ここが確定した時点で reader を固定する。

### Phase 1: 単一 Python スクリプト

本回答の `airgap_force_fft_prototype.py` を使用する。

目標:

- 実データが読める
- `sigma_r`, `sigma_t` が得られる
- torque, UMP が得られる
- Fr/Ftheta の 2D FFT が描ける

### Phase 2: Codex による整理

- モジュール分割
- dataclass / config 整理
- pytest 追加
- docstring / type hint
- エラー処理
- CLI 整理
- サンプルデータ追加

### Phase 3: FEM 結果による妥当性確認

最低限、

1. FEM のトルク波形
2. `sigma_t` 積分トルク
3. `Ftheta` の `m=0` 成分

を比較する。

また、

1. 直接積分 UMP
2. `Fr` の `m=1` 成分

の整合を確認する。

### Phase 4: GL80 評価

GL80 について、

- 無負荷
- 負荷
- ピークトルク近傍
- 複数ロータ位置

で評価する。

確認対象:

- torque ripple
- cogging torque
- `Fr(m=1)`
- `Fr(m=2)`
- その他の支配的低空間次数
- `Ftheta` 支配モード

### Phase 5: 偏心解析

理想対称モデルでは UMP がほぼゼロになる可能性が高い。

軸受・シャフトへの電磁横荷重を評価する場合は、

- static eccentricity
- dynamic eccentricity
- 磁石位置ずれ
- 磁化ばらつき

などを別ケースとして解析する。

特に `m=1` を評価する場合、
周期対称セクタモデルから全周を複製しただけでは
本来の非対称 UMP を生成できない点に注意する。

---

## 13. ロボット用途への拡張

サービスロボットでは一定回転速度運転だけを想定しない。

### 13.1 定常速度解析

通常の

```text
theta x time
```

2D FFT を使用できる。

### 13.2 変速運転

全時間を一括 FFT するとスペクトルが広がる。

将来は以下を追加する。

- rotor angle domain resampling
- order tracking
- STFT
- windowed 2D FFT

優先順位としては、
ロボット用途では rotor angle domain の order tracking を先に実装する。

### 13.3 停止・保持トルク

停止時は時間 FFT より、

- rotor position
- current (`id`, `iq`)
- spatial order

をパラメータにして

```text
Fr(m; rotor_angle, id, iq)
Ftheta(m; rotor_angle, id, iq)
```

を評価する方が有用。

---

## 14. 構造・振動解析との接続

将来、電磁力を構造解析へ渡す際には
単なる FFT 振幅だけでなく複素位相を保存する。

必要データ:

```text
m
frequency/order
complex amplitude
phase
```

これにより

- 進行波方向
- 構造モードとの対応
- 複素加振力
- モード重ね合わせ

へ発展できる。

電磁界解析側でギャップ上の位置情報を保持しておけば、
さらに節点力へのマッピングも可能。

---

## 15. 最適化との接続

将来的な目的関数候補:

```text
maximize average / peak torque
minimize cogging torque
minimize torque ripple
minimize Fr(m=1)  # UMP
minimize Fr(m=2)  # ovalizing force
minimize selected low-order radial force harmonics
```

構造解析まで接続した後は、

```text
minimize vibration acceleration
minimize housing displacement
minimize bearing electromagnetic load
```

へ拡張する。

---

## 16. Codex への初回指示案

以下を Codex に依頼する。

1. `airgap_force_fft_prototype.py` を基準実装として読む。
2. 数式・FFT 規約・出力値を変更せずモジュール化する。
3. synthetic test を先に実装する。
4. `m=0 Ftheta -> torque` を検証する。
5. `m=1 Fr -> Fx/Fy/UMP` を検証する。
6. 正負の spatial order / temporal frequency を保持する。
7. FFT 振幅定義を README に明記する。
8. 360° 全周版を最初に完成させる。
9. sector reconstruction や不等間隔補間は別 PR とする。
10. GL80 実データ reader は入力形式確定後に追加する。
