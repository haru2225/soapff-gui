# soapff-GUI：粘土鉱物の初期構造を選ぶだけで LAMMPS を回す

構造ファイル（元素・座標・箱だけ）を選ぶと、同梱の soapff-cui が SOAP で全原子の電荷と ClayFF の型を決め、LAMMPS の入力（`system.data`、`system.settings`、エネルギー最小化までの `in.lammps`）を作ります。その下に MD（fix、compute、run など）を一覧から選んで足し、「計算を開始」でこのコンピュータか SSH の先の LAMMPS を回し、解析（拡散係数、SOVA の PDF・RDF・配位数・S(Q)、compute の結果、debyer の XRD）と OVITO での表示までを行います。力場は ClayFF（Cygan et al. 2004 の表 1・表 2）と柔軟な SPC 水です（時間刻み 0.5 fs）。

SOAP のライブラリ（Na モンモリロナイトと 12 のバルクの結晶による SOAP と電荷の一対一の対応、13 構造、5,878 原子）は `src/soapff_gui/data/library_clayff_bulk.npz` に同梱していて、指定は不要です。別のライブラリは「計算設定」の「ライブラリと詳細」で選べます。

## 例

`examples/montmorillonite_2W/`：Na モンモリロナイト 2 水和層（2,502 原子）の初期構造と、GUI での手順、期待される結果（力場の作成と MD）、`diffusion/` に論文の拡散係数の計算（300、330、366 K、Greathouse ら 2016 との比較）。

`examples/library_nvt_npt/`：ライブラリにある 13 の鉱物（Na モンモリロナイトと 12 の結晶）それぞれの構造ファイルと、GUI で作った NVT と NPT の計算（入力、出力、解析）。

`examples/muscovite_2M1/`：白雲母 2M1（672 原子）の 0 K の弾性定数と 300 K の引張・圧縮・せん断（論文の計算）。物性の計算（＋ 物性）は今の版にはなく、この例は削除前の版（コミット 2c43d80）で作ったものです。

## 計算できる鉱物と構造

soapff-cui の README の「計算できる鉱物と構造」と同じです。要点：
- Na モンモリロナイト 2 水和層（JAEA の入力）と、次のバルクの結晶（水なし）：石英、アルバイト、ブルーサイト、パイロフィライト（2 つの出典）、白雲母（2M1、1M）、金雲母（1M）、nanpingite（Cs 雲母）、カオリナイト、リザーダイト、タルク。どれも、ライブラリと違う大きさの超格子（216〜720 原子）と、300 K の MD（NVT と NPT のあと）の構造で、全原子の型と電荷が正解（結合のつながりから決めた ClayFF の型）と一致しました。
- 力場を作るとき、SOAP で決めた型を、結合のつながりから決めた ClayFF の型と照らし合わせます（H のつき方が同じ型どうし。架橋 O、OH の O、陽イオン）。ライブラリで決まらなかった原子（fallback）か、型が食い違う原子（例：obts が ob になる）があると、使ったライブラリのパスとともに警告を出します（`report.json` の `type_check`、`charges.csv` の列 `bond_type`）。ライブラリにその鉱物の環境がないときに起こるので、計算設定 › ライブラリと詳細 の「同梱のものに戻す」で同梱のライブラリを選ぶか、参照を加えてください。
- annite（Fe²⁺ の雲母）は入れていません。ClayFF（Cygan et al. 2004 の表 1）には八面体の Fe³⁺ しかないためです。
- 金雲母は、ライブラリから自分の参照を抜くと、四面体の Al と 3 つの Mg に接する O（16 原子）が ob になります（ほかの鉱物に同じ環境がない）。新しい種類の鉱物は、その鉱物の参照をライブラリに入れてください。
- 表面、鉱物と水の界面（水酸基で終端した表面を含む）、エッジは、ライブラリにその環境がないので扱えません。
- 力場は ClayFF（原論文の表のとおり）と柔軟な SPC 水に固定です。

## インストール（GitHub から）

```
git clone git@github.com:haru2225/soapff-gui.git
sh soapff-gui/install.sh
```

これで、専用の Python の環境（`~/.soapff-gui/venv`）に soapff-gui と必要なライブラリ（numpy、scipy、ASE、dscribe、matplotlib）が入り、macOS では `~/Applications/soapff-GUI.app`（アイコンつき）ができます。Finder や Launchpad から「soapff-GUI」で起動できます（Dock にドラッグしておくと便利です）。

- Python 3.9 以上で tkinter（Tk 8.6 以上）のあるものが要ります。`python3`、`python3.13`〜`python3.9`、Homebrew、python.org の Python から自動で探します。Homebrew の Python には tkinter がないことが多く、そのときは uv（`brew install uv`）があれば uv の Python（tkinter つき）で環境を作ります。どちらもないときは `brew install python-tk@3.12` か `brew install uv` のあと、もう一度実行してください。使う Python は `PYTHON=/path/to/python sh soapff-gui/install.sh` でも選べます。
- LAMMPS がもう入っている Python の環境（conda や pixi など）にそのまま入れるときは、`PYTHON=その環境の python sh soapff-gui/install.sh --here` とします。その環境の `lmp` と `mpirun` を自動で使います。
- LAMMPS（`lmp`）と、並列にするなら `mpirun` は別に必要です（SSH で別のコンピュータを使うときは、そのコンピュータの LAMMPS を使います）。見つからないときは「計算設定」の「実行先」で場所を指定します。
- 解析の PDF・RDF・S(Q) には SOVA（sovapy）が、XRD には Debyer（https://github.com/wojdyr/debyer をビルドした `debyer`）が要ります。どちらも「計算設定 → 解析」で場所を指定できます（SOVA は `~/sova-cui` や `~/*/sova-cui`、debyer は PATH や `~/*/debyer/debyer/debyer` を自動で探します）。
- コマンドで使う `soapff-cui`（`soapff-cui make 構造ファイル -o 出力 --library …`）も一緒に入ります。力場の作成と解析の本体は `src/soapff_cui/` です。
- 更新は `cd soapff-gui && git pull` だけです（`-e` で入れているので、入れ直しは要りません）。

## 起動

```
soapff-gui              # = python -m soapff_gui
soapff-gui --install-app   # macOS：~/Applications/soapff-GUI.app を作る（Finder・Launchpad・Dock から起動できる）
```

`soapff-GUI.app` は、入れた Python で soapff-gui を起動するだけの小さなアプリです（AppleScript のアプレット。アイコンは `src/soapff_gui/data/icon.png`）。Python の環境を入れ替えたら、もう一度 `--install-app` を実行してください。

LAMMPS（`lmp`、`lmp_mpi`、`lmp_serial`）と `mpirun` は、Python と同じ環境の `bin` を先に、次に PATH を探します（`lmp` と同じフォルダの `mpirun` を優先）。見つからないときは「計算設定」の「実行先」で場所を指定してください。

## 画面

| 場所 | 内容 |
|---|---|
| メニュー | ファイル（新規計算、構造を変更、既存の出力フォルダを開く）、編集（コマンドを挿入、行を削除、上へ・下へ移動、元に戻す）、表示（解析結果、OVITO、計算履歴、構造ライブラリ）、計算（計算を開始、停止、入力を作り直す、解析、計算設定）、ヘルプ |
| ツールバー | 「新規計算」で構造ファイル（extxyz、cif、POSCAR、LAMMPS データ。元素・座標・箱だけで足ります）を選ぶと、soapff-cui が力場と入力（`in.lammps`）を作ります（出力フォルダは構造ファイルの隣の `<名前>_soapff/`）。「計算設定」で MD の条件・実行先・ライブラリを設定し、「実行・解析」で解析結果を表示します。右の「計算を開始」で LAMMPS を回します（計算中は「停止」） |
| 左の欄 | プロジェクト（新しい計算、計算履歴、構造ライブラリ）と、最近使った計算（クリックで開く） |
| 入力ファイル | 構造の名前と原子数、SOAP の割り当ての要約（粘土・水・イオンの数、SOAP / fallback の数、最小の類似度）、出力フォルダ。その下に `in.lammps`、`system.settings`（力場の設定）、`system.data`（原子・結合・角度。読み取り専用で、行を選ぶと Atoms なら原子番号・分子番号・型（名前つき）・電荷・座標のように列の名前つきで表示）をタブで切り替えて 1 行ずつ表示します。構造ファイルを入れた直後の `in.lammps` は**エネルギー最小化までの基本部分**だけです（soapff-cui の `make --minimize-only`。前の版で作った入力は「編集 → minimize より後をすべて消す」で同じ形にできます）（力場の読み込み、近接リスト、出力、最小化。周期でない面の壁）。MD はその下に、「＋ fix」「＋ compute」「＋ その他」から LAMMPS のコマンドを選んで足します。行をクリックして選び、上下キーで移動、Delete で削除。直した内容はすぐ保存されます（自動保存）。右上の「境界条件」で x・y・z を選ぶと（system.settings の `boundary` の行でも可）、力場と入力を作り直します |
| 選択したコマンド | 選んだ行の値を、名前付きの欄（fix nvt なら ID、対象グループ、スタイル、開始温度、終了温度、減衰）で直せます。下に「生成される行」を表示します |
| 実行 | 実行先（このコンピュータ、または登録したコンピュータ）、「計算を開始」、状態（段階と進み具合）、温度・圧力・ポテンシャルエネルギー、実行ログ。「解析」「OVITO で開く」「フォルダ」 |
| 状態の帯 | 構造と原子数、実行先 |

進み具合は `in.lammps` の `run` の合計ステップに対する割合です（行を直すと自動で変わります）。

### fix・compute の追加

「＋ fix」「＋ compute」「＋ その他」で、LAMMPS のコマンドの一覧から選んで、選んだ行のあとに入れます。名前や説明で絞り込め、入る行と欄の名前を確かめてから「挿入」します。「LAMMPS のマニュアル」で、そのコマンドの説明のページ（docs.lammps.org）を開きます。ID は、同じ ID があれば自動で番号を付けます。温度・圧力・時間刻み・ステップ数は、計算設定の MD の条件と in.lammps の timestep から入ります。

| 一覧 | 入っているもの |
|---|---|
| fix | 時間積分・温度・圧力（nvt、npt（等方、z、aniso、tri）、nphug（衝撃圧縮のユゴニオ状態）、nph、nve、langevin、temp/csvr、temp/berendsen、temp/rescale、press/berendsen、box/relax、rigid）、拘束・束縛（shake、rattle、setforce、spring/self、momentum、recenter、restrain、viscous）、外場・壁（efield、addforce、aveforce、wall/reflect、wall/lj93、deform、heat）、出力・平均（ave/time、ave/chunk、ave/correlate、ave/atom、dt/reset） |
| compute | 拡散・構造（msd、msd/chunk、vacf、rdf、displace/atom、coord/atom、cluster/atom、centro/atom、orientorder/atom、gyration、com、chunk/atom）、温度・エネルギー・圧力（temp、temp/com、temp/partial、pe、pe/atom、ke、ke/atom、pressure、stress/atom、heat/flux、group/group、pair）、原子・結合の量（property/atom、reduce、count/type、bond、angle、bond/local、angle/local、dipole） |
| その他 | run、timestep、velocity、unfix、uncompute、reset_timestep、minimize、thermo、thermo_style、dump（OVITO 用）、write_data、write_restart、group |
| soapff の手順 | 拡散係数の MSD（「解析」で使う msd_<名前>.dat を書く）、soapff の緩和（初速度、時間刻み 0.5 fs、緩和 NVT、ステップのリセットまで。本計算は一覧から足す） |

`in.lammps` に run がないまま「計算を開始」を押すと、最小化だけでよいか確認します。

### 計算設定

| タブ | 内容 |
|---|---|
| MD の条件 | 温度、緩和 NVT（最小化のあと、記録しない）、NVT、NPT、圧力、出力間隔、圧力の方向、境界条件。次に入力を作るとき（新規計算、「保存して入力を作り直す」）に使います |
| 実行先 | このコンピュータの lmp、mpirun、並列数。登録したコンピュータ（SSH）の一覧と、追加・編集・削除・接続テスト |
| ライブラリと詳細 | SOAP ライブラリ、類似度のしきい値、カットオフ |
| 解析 | 終わったら解析するか、行う解析（拡散係数、SOVA、compute の結果、XRD）、平衡化として除く時間、XRD の条件、SOVA と debyer の場所（下の「解析」） |

## 別のコンピュータで計算する（SSH）

「計算設定」→「実行先」で、コンピュータを名前を付けて登録し、「実行」の欄の実行先でその名前を選んで「計算を開始」を押します。

1. 登録する項目：名前（例：研究室のサーバー）、ホスト、ユーザー、ポート、相手側のフォルダ（計算を置く場所。既定 `~/soapff_runs`）、相手側の `lmp` と `mpirun`（名前かフルパス）、最初に実行するコマンド（例：`module load lammps`、`source ~/.bashrc`。SSH で実行するコマンドは `~/.bashrc` を読まないことが多いため）、並列数。
2. **パスワードは保存しません。** 計算を始めるたびと、接続テストのたびに聞きます。空欄のまま OK にすると、鍵でログインします。
3. 計算を始めると、このコンピュータで力場と入力を作り、`system.data`、`system.settings`、`in.lammps` を相手側の `<フォルダ>/<出力フォルダの名前>/` に送って、LAMMPS を `nohup` で起動します（接続が切れても計算は続きます）。画面には相手側の画面出力を流し、進み具合を表示します。
4. 終わると、相手側の計算フォルダを中身ごと（`scp -r`、`cp -r` と同じ）出力フォルダに受け取り、解析はこのコンピュータで行います。
   「ファイル → 別のコンピュータからフォルダを受け取る」で、相手側の好きなフォルダ（ほかの計算や、相手側で作ったデータ）を一覧から選んで、丸ごとこのコンピュータにコピーできます。in.lammps のあるフォルダは「既存の出力フォルダを開く」で開けます。
5. 「停止」は相手側の mpirun と LAMMPS を止めます。計算が終わるか止めると、接続を閉じます（次はまたパスワードを聞きます）。

使うのはこのコンピュータの `ssh` と `scp` だけです。相手側で使うのは、相手のコンピュータにもともと入っている LAMMPS（`lmp`）と `mpirun` だけで、SOAP のライブラリと力場の作成（soapff-cui）、解析（SOVA）、OVITO はこのコンピュータで動きます（相手側に Python、soapff、OVITO は要りません）。相手の `lmp` には、送った入力の機能（`lj/cut/coul/long`、`pppm`、`harmonic`。周期でない方向があれば `msm`）が要ります。
並列にするときは、`lmp` と同じ MPI で作られた `mpirun` が必要です（違うと、1 並列の LAMMPS が並列数の分だけ別々に動きます）。接続のたびに `lmp -h` の MPI と `mpirun` の種類を調べ、`lmp` と同じフォルダに `mpirun` があればそれを使い、合わないときは警告します。相手側は Linux（または macOS）の sh を想定しています。初めて接続するホストの鍵は自動で `~/.ssh/known_hosts` に登録します。

## 解析

計算が終わると（「計算設定 → 解析」で「計算が終わったら解析する」のとき）、または「解析」ボタンで、選んだ解析を行います。「解析」ボタンでは、行う解析と条件を選んでから始めます。SOVA と XRD に使う軌道は、出力フォルダの `*.lammpstrj`（既定は MD の `traj.lammpstrj`）から選べます。`traj.lammpstrj` がなくても、MSD や fix ave/* の出力があれば解析できます。結果は「解析結果」のウィンドウに、まとめ（数値）と図のタブで表示します。解析はこのコンピュータで行います（SSH で計算したときは、受け取ったフォルダで）。

| 解析 | 使うもの | 内容 | ファイル（`analysis/`） |
|---|---|---|---|
| 拡散係数 | LAMMPS の compute msd（カタログの「拡散係数の MSD」が書く `msd_<名前>.dat`） | アインシュタインの式 D = 傾き / (2d)。平衡化のあとを直線近似。全体・面内 xy・垂直 z。有限サイズ補正なし | `msd_fit.csv`、`analysis.png`（左） |
| SOVA | SOVA（sovapy）：atom_pair_hist、pair_dist_func、radial_dist_fun、partial_structure_factor、xcoeff、structure_factor_xray | 水O–水O、イオン–水O、イオン–鉱物O と全元素の組の PDF g(r)、RDF N(r) = 4πr²ρg(r)、第 1 ピーク、配位数（N(r) を第 1 極小まで積分）。全原子の X 線の構造因子 S(Q)（g(r) を箱の幅の半分まで使う。Q < 2π/r_max は決まらないので図で網掛け。層状・界面の系では g(r) が 1 に近づかないので、低 Q の値は使わない） | `sova_pdf.csv`、`sova_rdf.csv`、`sova_sq.csv`、`sova_summary.json`、`analysis.png`（中・右）、`sova_sq.png` |
| compute の結果 | `in.lammps` の `fix ave/time`、`ave/chunk`、`ave/correlate`、`ave/histo` のうち `file` で書き出したもの（カタログの compute を fix ave/* で書き出したもの） | 1 行 1 時刻のファイル：列ごとに時間変化の図と、平衡化のあとの平均。ブロックのファイル（mode vector、chunk、ヒストグラム、相関）：平衡化のあとの出力を平均した図（列ごとに別の段） | `compute_<ファイル名>.png`、`.csv`、`_average.csv` |
| XRD | Debyer（デバイの散乱式）。`debyer -x -l λ` | 最後のフレーム（既定 3 つ）の粉末 X 線回折パターン I(2θ) の平均と、強いピークの 2θ・d・相対強度。**結晶子**（既定。結晶・粘土鉱物）：各辺が「結晶子の大きさ」（既定 50 Å）以上になるようにセルを並べた有限の結晶で計算（周期境界なし）。ピークは結晶子の大きさで広がり、層の少ない方向の低角のピーク（001 など）は低角側にずれます。d が結晶子の 1/3 を超える極大は大きさによる振動として別に示します。**周期境界**（液体・ガラス。直方体の箱のみ）：カットオフ = 箱の短い辺の半分、sinc 減衰、箱の数密度。計算時間は原子数の 2 乗に比例します（2 万原子で 1 フレーム約 40 秒） | `xrd.csv`、`xrd.png` |

数値のまとめは `summary.txt`（`summary.json`）です。コマンドでは `soapff-cui analyze 出力フォルダ [--what msd,sova,compute,xrd] [--traj traj.lammpstrj] [--skip-ps 10] [--debyer パス] [--wavelength 1.5406] [--two-theta 3,80,0.02] [--xrd-frames 3] [--xrd-mode cluster|periodic] [--xrd-size 50] [--sova-python パス] [--sova-home パス]` です。

## 出力（出力フォルダ）

soapff-cui の出力（`system.data`、`system.settings`、`in.lammps`、`charges.csv`、`report.json`）に加えて、LAMMPS の `log.lammps`、`traj.lammpstrj`（出力間隔ごと、既定 1 ps）、`final.data`（最後の構造）。MSD の `msd_<group>.dat`（LAMMPS）。解析を行うと `analysis/` に結果ができます（下の「解析」）。

## 構成

```
src/soapff_gui/   app.py           ウィンドウ（tkinter）。計算は別スレッドで行い、画面はキューで更新
                  lammps_input.py  in.lammps を行ごとに扱う（分解、欄の名前、最小化までの基本部分、ID、進み具合）
                  catalog.py       追加できる LAMMPS の fix・compute・その他のコマンドの一覧
                  remote.py        SSH で別のコンピュータを使う（ssh / scp、nohup での起動と停止、受け取り）
                  runner.py        同梱の soapff-cui、LAMMPS・mpirun・SOVA の場所、LAMMPS の起動
                  data/            SOAP のライブラリ（library_clayff_bulk.npz）
src/soapff_cui/   同梱の soapff-cui（構造ファイル → ClayFF の LAMMPS 入力、解析）。_clayff は型付けの部品、
                  analysis.py は解析（MSD とイオン伝導度、fix ave/* のファイル、debyer）、sova_structure.py は SOVA の部分
examples/         例：ライブラリの 13 の鉱物の NVT と NPT、Na モンモリロナイト 2 水和層（力場、MD、拡散係数）、
                  白雲母 2M1（論文の弾性定数、引張、せん断。削除前の版で作ったもの）
```

## 限界

- 力場とライブラリの限界は soapff-cui と同じです（上の「計算できる鉱物と界面構造」）。
- このコンピュータか、SSH で入れるコンピュータで、LAMMPS を直接起動します。PBS などのジョブ投入はしません。
