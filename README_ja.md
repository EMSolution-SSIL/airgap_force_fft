# Air-Gap Force FFT

モータのエアギャップ磁束密度結果を後処理する Python ツールです。`Br` と `Btheta` を読み込み、Maxwell 応力を計算し、積分力およびトルクを評価するとともに、2次元 FFT により時間・空間高調波モードを抽出します。

同梱のサンプルは、GL80 PMSM の健全位置（healthy position）における解析結果です。

本プロジェクトは [MIT License](LICENSE) の下で公開されているオープンソースソフトウェアです。

## リポジトリ構成

```text
.
|-- src/                         # 処理スクリプトおよび再利用可能なモジュール
|   |-- airgap_force_fft.py       # CLI エントリーポイント
|   |-- airgap_io.py              # JSON/CSV リーダーおよびセクタ展開
|   |-- airgap_fft.py             # Maxwell 応力、FFT、トルク比較
|   `-- airgap_plotting.py        # プロット用ユーティリティ
|-- tests/                       # 単体テスト
|-- data/
|   `-- GL80/healthy_position/   # 公開サンプル入力および生成結果
|-- docs/                        # 検討メモおよび解析レポート
|-- config.yaml                  # GL80 サンプルのデフォルト設定
|-- pyproject.toml               # パッケージ設定
`-- requirements.txt             # 実行時依存パッケージ
```

## クイックスタート

PyPI から公開パッケージをインストールします。

```powershell
python -m pip install airgap-force-fft
```

開発用途では、仮想環境を作成してリポジトリの依存パッケージをインストールします。

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

デフォルトの GL80 サンプルを実行します。

```powershell
python .\src\airgap_force_fft.py --config .\config.yaml
```

主な出力は、デフォルトでは `data/GL80/healthy_position/fft_output` 以下に保存されます。

## 入力形式

以下のエアギャップ磁束密度入力に対応しています。

- EMSolution/pyemsol の `gapB` JSON
- `time`、`theta_deg`、`Br`、`Btheta` などの列を持つ long-form CSV

EMSolution JSON では、ファイルに含まれている場合、形状情報および対称性メタデータを自動的に読み取ることができます。CSV 入力では、不足するメタデータを `config.yaml` に明示的に設定してください。特に以下の項目が必要です。

- `input.sector.periods`
- `input.sector.symmetricity`
- `geometry.radius_m`
- `geometry.stack_length_m`

360度全周の CSV データを使用する場合は、`input.sector.periods: 1` および `input.sector.symmetricity: 0` を指定してください。

## トルク参照ファイル

トルク比較では、以下の2種類の参照形式に対応しています。

- 従来の `transient_results.json` 形式（`Time.data` および `Torque.*.data`）
- pyemsol の `output_transient.json`。`postData.forceNodal.forceNodalData` のうち、`propertyNum` が `stator` または `rotor` のエントリから `forceMZ` を読み込みます。

比較プロットでは、ロータ正方向の符号規約として `-Ftheta integrated`、`-stator`、`rotor` を使用します。

## レポート

主要な解析レポートは `docs/` にあります。

- [GL80 ベースライン解析](docs/GL80_analysis_report.md)
- [初期検討メモ](docs/airgap_force_fft_discussion.md)
- [実装計画](docs/airgap_force_fft_implementation_plan.md)

## テスト

単体テストは以下のコマンドで実行します。

```powershell
python -m unittest discover -s tests -v
```

## パッケージング

開発用 editable install:

```powershell
pip install -e .
```

インストール後は、CLI を以下のように実行することもできます。

```powershell
airgap-force-fft --config config.yaml
```

## 英語版README

英語版は [`README.md`](README.md) を参照してください。
