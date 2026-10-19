# soapff-GUI：粘土鉱物の初期構造を選ぶだけで LAMMPS を回す

[soapff-cui](../soapff-cui) のウィンドウ版です。初期構造のファイルを選んで「実行」を押すと、

1. soapff-cui が出力フォルダに LAMMPS の入力一式（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）を書き出し、
2. そのまま LAMMPS（`mpirun -np N lmp -in in.lammps`）を起動します。LAMMPS は、最小化 → 10 ps の緩和 NVT（記録しない）→ 本計算の NVT（→ 指定すれば NPT）の順に、ClayFF（Cygan et al. 2004 の表 1・表 2）と柔軟な SPC 水、時間刻み 0.5 fs で回ります。

入力ファイルを開いたり書き換えたりする必要はありません。SOAP のライブラリは、隣の soapff-cui の `examples/minerals/library_clayff_bulk.npz`（Na モンモリロナイトと 4 つの鉱物のバルクの結晶による SOAP と電荷の一対一の対応、5 構造、3,432 原子）を自動で使うので、指定は不要です。それが見つからないときは、GUI に同梱したコピー（`soapff_gui/data/library_clayff_bulk.npz`）を使います。以前の設定で別のライブラリを保存しているときは、「詳細設定」で変えてください。

## 例

`examples/montmorillonite_2W/`：Na モンモリロナイト 2 水和層（2,502 原子）の初期構造と、GUI での手順、期待される結果（力場の作成と MD）。

## 計算できる鉱物と構造

soapff-cui の README の「計算できる鉱物と構造」と同じです。要点：
- Na モンモリロナイト 2 水和層（JAEA の入力）と、石英・アルバイト・ブルーサイト・パイロフィライトのバルクの結晶（水なし）。ライブラリより大きい超格子の結晶（216〜720 原子）で全原子の型と電荷が正解と一致しました。
- 表面、鉱物と水の界面（水酸基で終端した表面を含む）、エッジは、ライブラリにその環境がないので扱えません。
- 力場は ClayFF（原論文の表のとおり）と柔軟な SPC 水に固定です。

## インストール（GitHub から）

```
# 1. soapff-cui（力場を作る本体）
git clone git@github.com:haru2225/soapff-cui.git && pip install -e soapff-cui
# 2. soapff-gui
git clone git@github.com:haru2225/soapff-gui.git && pip install -e soapff-gui
soapff-gui              # = python -m soapff_gui
```

Python には tkinter が必要です（多くの配布に含まれます）。LAMMPS（`lmp`）と `mpirun` は別に必要です。

## 起動

```
# A. macOS：Finder で soapff-gui.command をダブルクリック
#    （Python は既定で ~/ScoreMD/.pixi/envs/default/bin/python。PY で変えられます）

# B. ターミナルから
PYTHONPATH=~/clay/soapff-gui ~/ScoreMD/.pixi/envs/default/bin/python -m soapff_gui
```

soapff-cui がインストールされていないときは、`$SOAPFF_HOME`、このフォルダの隣の `soapff-cui/`（`~/clay/soapff-cui`）、`soapff/`（古い名前）の順に探します。見つからないときは、起動したときに知らせます。
LAMMPS（`lmp`、`lmp_mpi`、`lmp_serial`）と `mpirun` は、Python と同じ環境の `bin` を先に、次に PATH を探します。見つからないときは「詳細設定」で場所を指定してください。

## 画面

| 項目 | 内容 |
|---|---|
| 初期構造 | extxyz（Lattice 付き）、cif、POSCAR、LAMMPS データ。元素・座標・箱だけで足ります。選ぶと原子数・元素・箱の大きさを表示します |
| 出力フォルダ | 既定は構造ファイルの隣の `<名前>_soapff/`。前の結果があるときは上書きを確認します |
| MD の条件 | 温度、緩和 NVT の長さ（既定 10 ps、0 でなし。最小化のあと、本計算の前に回し、記録しない）、本計算の NVT の長さ、NPT の長さ（0 でなし）、圧力、出力間隔（ps。拡散係数には 0.1〜0.5）、圧力の方向（界面・層間は `z`、バルク結晶は `aniso` / `tri`）、並列数（MPI） |
| 境界条件 | x・y・z それぞれ `p`（周期）、`f`（固定）、`m`、`s`、面ごとに違う `fm` など（LAMMPS の boundary と同じ）。既定は `p p p`。周期でない方向があるときは NPT は使えません。長距離クーロン力の計算法は soapff-cui が選びます（soapff-cui の README の「境界条件」） |
| 詳細設定 | SOAP ライブラリ、lmp と mpirun の場所、類似度のしきい値（既定 0.95）、カットオフ（既定 10 Å）、解析で除く平衡化の時間（空欄で最初の 20 %）、SOVA の Python とフォルダ（この Python に sovapy があれば空欄。なければ `~/sova-cui` や `~/*/sova-cui` を自動で探す） |
| ▶ 実行 / ■ 停止 | 力場の作成 → 最小化 → 緩和 NVT → 本計算の MD。停止は LAMMPS を止めます |
| 入力だけ作る | LAMMPS を回さず、入力一式だけ書き出します（スパコンに持っていくときなど） |
| 解析 / 終わったら解析 | 拡散係数（LAMMPS の `compute msd`）と、PDF g(r)・RDF（SOVA）。「終わったら解析」にチェックすると MD のあと自動で行います。「解析」ボタンは計算中でも、それまでの結果で解析します。結果のウィンドウにまとめと 3 枚組の図（MSD、PDF、RDF）を表示し、CSV のフォルダを開けます |
| OVITO で開く | 出力フォルダの `traj.lammpstrj`（なければ `final.data`、`system.data`）を OVITO で開きます。計算中でも、それまでのフレームを開けます。MSD と PDF・RDF 以外の解析（密度分布など）は OVITO の中で行えます。macOS は Ovito.app、それ以外は `ovito` コマンドが必要です |
| 状態 | 段階（力場作成・最小化・緩和 NVT・MD）、緩和と本計算を合わせた進み具合、ステップ・温度・圧力・ポテンシャルエネルギー・密度 |

力場を作ったあと、原子数、SOAP で電荷を決めた原子と fallback の原子の数、最小の類似度を表示します。fallback の原子や重なりの警告があるときは橙色で知らせます（ライブラリにない環境が含まれている可能性があります。`charges.csv` を確認してください）。

設定は `~/.config/soapff-gui/settings.json` に保存し、次回の起動で使います。

## 出力（出力フォルダ）

soapff-cui の出力（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）に加えて、LAMMPS の `log.lammps`、`traj.lammpstrj`（出力間隔ごと、既定 1 ps）、`final.data`（最後の構造）。MSD の `msd_<group>.dat`（LAMMPS）。解析を行うと `analysis/` に `summary.txt`、`summary.json`、`analysis.png`、`msd_fit.csv`（拡散係数）、`sova_pdf.csv`、`sova_rdf.csv`、`sova_summary.json`（SOVA）ができます（内容は soapff-cui の README の「解析」）。

## 構成

- `soapff_gui/runner.py`：GUI を使わない部分（soapff と LAMMPS を探す、soapff-cui の `make` を呼ぶ、LAMMPS を起動する、画面出力から段階と thermo の値と進み具合を読む）。
- `soapff_gui/app.py`：tkinter のウィンドウ。計算は別スレッドで行い、画面はキューで更新します。

## 限界

- 力場とライブラリの限界は soapff-cui と同じです（上の「計算できる鉱物と界面構造」）。
- 1 台の計算機で回すためのものです。PBS などのジョブ投入はしません（「入力だけ作る」で入力を作って持っていってください）。
