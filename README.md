# soapff-GUI：粘土鉱物の初期構造を選ぶだけで LAMMPS を回す

[soapff-cui](../soapff-cui) のウィンドウ版です。初期構造のファイルを選んで「実行」を押すと、

1. soapff-cui が出力フォルダに LAMMPS の入力一式（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）を書き出し、
2. そのまま LAMMPS（`mpirun -np N lmp -in in.lammps`）を起動します。

入力ファイルを開いたり書き換えたりする必要はありません。SOAP のライブラリは、隣の soapff-cui の `examples/minerals/library_clayff.npz`（いちばん新しいもの）を自動で使うので、指定は不要です。それが見つからないときは、GUI に同梱したコピー（`soapff_gui/data/library_clayff.npz`）を使います。

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
| MD の条件 | 温度、NVT の長さ、NPT の長さ（0 でなし）、圧力、圧力の方向（界面・層間は `z`、バルク結晶は `aniso` / `tri`）、並列数（MPI） |
| 境界条件 | x・y・z それぞれ `p`（周期）、`f`（固定）、`m`、`s`、面ごとに違う `fm` など（LAMMPS の boundary と同じ）。既定は `p p p`。周期でない方向があるときは NPT は使えません。長距離クーロン力の計算法は soapff-cui が選びます（soapff-cui の README の「境界条件」） |
| 詳細設定 | SOAP ライブラリ、lmp と mpirun の場所、類似度のしきい値（既定 0.95）、カットオフ（既定 10 Å） |
| ▶ 実行 / ■ 停止 | 力場の作成 → 最小化 → MD。停止は LAMMPS を止めます |
| 入力だけ作る | LAMMPS を回さず、入力一式だけ書き出します（スパコンに持っていくときなど） |
| OVITO で開く | 出力フォルダの `traj.lammpstrj`（なければ `final.data`、`system.data`）を OVITO で開きます。計算中でも、それまでのフレームを開けます。RDF、配位数、密度分布などの解析は OVITO の中で行います。macOS は Ovito.app、それ以外は `ovito` コマンドが必要です |
| 状態 | 段階（力場作成・最小化・MD）、MD の進み具合、ステップ・温度・圧力・ポテンシャルエネルギー・密度 |

力場を作ったあと、原子数、SOAP で電荷を決めた原子と fallback の原子の数、最小の類似度を表示します。fallback の原子や重なりの警告があるときは橙色で知らせます（ライブラリにない環境が含まれている可能性があります。`charges.csv` を確認してください）。

設定は `~/.config/soapff-gui/settings.json` に保存し、次回の起動で使います。

## 出力（出力フォルダ）

soapff-cui の出力（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）に加えて、LAMMPS の `log.lammps`、`traj.lammpstrj`（1 ps ごと）、`final.data`（最後の構造）。

## 構成

- `soapff_gui/runner.py`：GUI を使わない部分（soapff と LAMMPS を探す、soapff-cui の `make` を呼ぶ、LAMMPS を起動する、画面出力から段階と thermo の値と進み具合を読む）。
- `soapff_gui/app.py`：tkinter のウィンドウ。計算は別スレッドで行い、画面はキューで更新します。

## 限界

- 力場とライブラリの限界は soapff-cui と同じです（エッジ、K や Cs、ライブラリにない鉱物は正しく扱えません）。
- 1 台の計算機で回すためのものです。PBS などのジョブ投入はしません（「入力だけ作る」で入力を作って持っていってください）。
