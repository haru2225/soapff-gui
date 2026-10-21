# soapff-GUI

粘土鉱物の構造ファイル（元素・座標・箱だけ）を選ぶと、SOAP で全原子の ClayFF の型と電荷を決めて LAMMPS の入力を作り、MD を回して解析する GUI です。力場は ClayFF（Cygan et al. 2004 の表のとおり）と柔軟な SPC 水です。

## できること

- **力場の自動生成**：同梱のライブラリ（13 構造、5,833 原子）でいちばん似た原子の型と電荷を付け、結合・角度・LJ をそろえた `system.data`、`system.settings`、最小化までの `in.lammps` を書く。ライブラリにない環境の原子（類似度 0.95 未満）があると、力場を作らずに止まる。
- **MD の追加**：最小化より後は、LAMMPS の fix・compute・その他のコマンドの一覧から選んで加える（温度・圧力・時間刻み・ID は自動で埋める）。
- **スーパーセル**：OVITO の Replicate で a・b・c 方向に並べてから力場を作る。
- **膨潤曲線の計算**：層間に水分子を加えながら NPT を回し、各状態の箱を `swelling.dat` に記録する（Cygan et al. 2004 と同じ手順。図にするのは `examples/montmorillonite_swelling/analyze_swelling.py`）。
- **実行**：このコンピュータか、SSH で登録した別のコンピュータ（`ユーザー@ホスト` を入れるだけで登録）で LAMMPS を回す。
- **解析**：拡散係数（`compute msd` があるとき）、SOVA の PDF・RDF・配位数・S(Q)、Debyer の XRD。SOVA の dr、r の最大、Q の範囲と刻み、X 線／中性子、フレーム数、求める組と、XRD の波長・2θ・結晶子の大きさは「解析」の画面（計算設定 › 解析）で選べます。OVITO で表示。

## インストール

| OS | 手順 |
|---|---|
| macOS | `git clone https://github.com/haru2225/soapff-gui.git && sh soapff-gui/install.sh` → `~/Applications/soapff-GUI.app` から起動 |
| Linux | `sudo apt install python3-venv python3-tk`（Ubuntu の例）のあと、macOS と同じ。起動は `~/.soapff-gui/venv/bin/soapff-gui` |
| Windows | python.org の Python 3.10–3.12 を入れ、`git clone …` のあと `powershell -ExecutionPolicy Bypass -File soapff-gui\install.ps1` → デスクトップのショートカットから起動（動作は未確認。SSH の実行は使えないので、SSH を使うときは WSL2 に Linux の手順で入れる） |

- LAMMPS（KSPACE と MOLECULE パッケージ入り）と `mpirun` は別に入れます（例：`conda install -c conda-forge lammps openmpi`、`brew install lammps`）。見つからないときは「計算設定 › 実行先」で場所を指定します。
- **SOVA（sovapy、PDF・RDF・S(Q)）と Debyer（XRD）は `install.sh` が入れます**（SOVA は `pip install sovapy`（Python 3.11 以上）、Debyer は git から `~/.soapff-gui/bin/debyer` にビルド。macOS は Homebrew で `autoconf automake gengetopt libtool` を入れ、Linux はビルドツールが足りなければ入れ方を表示するので、それを実行してもう一度 `install.sh`）。Windows は SOVA だけ入り、Debyer は WSL2 で入れます。表示（OVITO、https://www.ovito.org）は別に入れます。
- Tk のある Python が見つからないときは `brew install python-tk@3.12` か `brew install uv` のあと、もう一度 `install.sh` を実行します。更新は `git pull` です（依存が増えたら `~/.soapff-gui/venv/bin/python -m pip install -e soapff-gui`）。

## 使い方

1. 「新規計算」（またはスーパーセルにするなら「スーパーセル」）で構造ファイルを選ぶ → 力場と最小化までの入力ができる。
2. 「＋ fix」「＋ compute」「＋ その他」から MD を加える（例：soapff の緩和 → NVT → NPT → 最後の構造を書き出す）。膨潤曲線は「＋ その他 › soapff の手順 › 膨潤曲線」。
3. 実行先を選んで「計算を開始」。終わると解析し、図と数値を表示する（`analysis/summary.txt`）。

## 別のコンピュータ（SSH）で計算する

力場の作成と解析はこのコンピュータで、LAMMPS だけを相手側（Linux か macOS）で動かします。

1. **相手側**：SSH でログインでき、`lmp`（KSPACE と MOLECULE 入り）と `mpirun` があること。
2. **登録**：「実行」の欄の「実行先」の横の **＋登録** を押し、`ssh` の接続先（例：`taro@192.168.0.10`、ポートが 22 でなければ `taro@192.168.0.10:10022`）を入れて「登録」。パスワードを聞かれるので答えると、接続して相手の CPU と LAMMPS を調べ、並列数（物理コア数）と `mpirun` を自動で設定します。
3. **実行**：実行先でそのコンピュータを選んで「計算を開始」。入力を送って `nohup` で起動し（接続が切れても続く）、終わると結果のフォルダを丸ごと受け取ります。
4. **削除**：実行先で選んで **－削除**。細かい設定（相手側のフォルダ、`lmp` と `mpirun` の場所、`module load lammps` など最初に実行するコマンド、並列数）は「計算設定 › 実行先」で変えられます。

パスワードは保存しません（登録、計算、接続テストのたびに聞きます）。鍵（`ssh-copy-id`）でログインできるなら、パスワードは空欄で OK です。

## 例（`examples/`）

| フォルダ | 中身 |
|---|---|
| `library_nvt_npt/` | ライブラリの 13 構造それぞれの NVT と NPT（GUI で作った入力・出力・解析） |
| `montmorillonite_2W/` | Na モンモリロナイト 2 水和層の拡散係数（300・330・366 K、Greathouse et al. 2016 と比較） |
| `montmorillonite_swelling/` | Cygan et al. 2004 と同じ Na モンモリロナイトの構造の作り方と膨潤曲線 |
| `muscovite_2M1/` | 白雲母 2M1 の弾性定数と引張・せん断（論文の計算） |

## 扱える構造と限界

- ライブラリ：Na モンモリロナイト（2 水和層、乾燥）、石英、ブルーサイト、パイロフィライト、白雲母（2M1、1M）、金雲母、nanpingite、カオリナイト、リザーダイト、タルク。これらの結晶は、大きさの違う超格子でも、300 K の MD のあとの構造でも、全原子の型と電荷が正解と一致しました。
- 表面、界面、エッジ、ほかの鉱物は、ライブラリにその環境がないので扱えません（その参照をライブラリに加える）。
- 力場は ClayFF と柔軟な SPC 水に固定です。

## 構成

```
src/soapff_gui/  GUI（app.py）、コマンドの一覧（catalog.py）、SSH（remote.py）、実行（runner.py）、ライブラリ（data/）
src/soapff_cui/  力場の作成（core.py）、解析（analysis.py）。コマンドの soapff-cui
examples/        上の例
```
