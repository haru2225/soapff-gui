# 検証の記録（論文 `paper/main.tex` の結果）

論文の数値は、ここにあるスクリプトと結果から出しています。軌道（`*.lammpstrj`）、`read_dump` 用のダンプ（`*.dump`）、検証用に作った一時ライブラリ（`*.npz`）は大きいので入れていません（スクリプトで作り直せます）。スクリプトのパスは作成者の計算機（`/Users/harutokono/clay/...`、元は `soapff-cui/examples/`）のままです。

| フォルダ | 内容 | 論文 |
|---|---|---|
| `library/` | ライブラリの参照：バルクの結晶の作成（`make_bulk.py`、`make_bulk_sheet.py`。層状ケイ酸塩 8 つは共著者提供の POSCAR から、元素の行のない 3 つは配位数で元素を確認）、ClayFF の型付け（`build_clayff_refs*.py`）、参照のデータ（`refs_clayff/`、13 構造）、ライブラリの中の一対一の確認（`library_clayff_sheet_check.json`、5,878 / 5,878）。ライブラリそのものは `src/soapff_gui/data/library_clayff_bulk.npz` | 「ライブラリの参照」、表「ライブラリの中での一対一の確認」 |
| `earlier_tests/` | 5 構造のライブラリで行った試験（JAEA の初期構造、4 つの結晶の超格子、界面）のスクリプト | 「構造ファイルから作った力場」「ライブラリより大きい結晶」 |
| `sheet_silicates/` | 層状ケイ酸塩 8 つ：違う大きさの超格子と 300 K の MD（NVT と NPT）のあとのフレームの型付け（`md_sheet.py`、`analyze_sheet.py` → `results_sheet.json`）、1 つ抜いたライブラリ（`loo_sheet.py` → `results_loo.json`）、前の試験の 13 構造のライブラリでのやり直し（`run_old_tests.py` → `results_old_tests.json`）、MD の入力と log | 「層状ケイ酸塩の結晶」、表「層状ケイ酸塩の結晶での試験」と MD の表 |
| `muscovite/` | 白雲母 2M1（672 原子）の力学物性。`muscovite_2M1_mech_gui/`：soapff-gui で構造ファイル → 「＋ 物性」（0 K の弾性定数、x y z の引張、xy xz yz のせん断）→ SSH の計算機（8 並列）で回し、GUI が解析したもの（`analysis/summary.txt`、図）。`e0test/`：0 K の弾性定数のひずみ（0.005 と 0.01、2 並列と 4 並列）の確認。`run/`、`tension_z16/`：ローカル（4 並列）で同じ入力を回したもの（再現性）。`muscovite_2M1_221_soapff/`：古い 5 構造のライブラリで作ってしまった力場での計算（obts が ob になった例。型の照合を入れるきっかけ） | 「白雲母の弾性定数と力学試験」（方法と結果）、白雲母の表と図 |
| `diffusion_T2/` | Na モンモリロナイト 2 水和層の拡散係数（300、330、366 K、4 区間ずつ）の入力、log、MSD、`analyze_diffusion.py` → `diffusion_summary.json` | 「拡散係数の先行研究との比較」と図 |
| `code_checks/argon/` | 有限温度の弾性定数（LAMMPS examples/ELASTIC_T/BORN_MATRIX/Argon に本ツールの出力を加えた `in.test`、`log.test`）と、グリーン–久保の粘度・熱伝導率の単位（`in.gk`、`log.gk`）を LAMMPS 自身の値と比べたもの | 「白雲母の弾性定数と力学試験」の「計算コードの確認」 |
| `code_checks/quartz/` | 石英の 0 K の弾性定数と、箱を等方に伸縮した圧力・エネルギーからの体積弾性率（`bcheck.in`、`bcheck.log`） | 同上 |
| `code_checks/type_check.py` | SOAP の型と結合のつながりの型の照合（誤検出 0、obts → ob の検出） | 「力場の自動生成」の割り当て |

主な結果：
- ライブラリ（13 構造、5,878 原子）の中で一対一の対応が成り立つ（l_max = 3、4、5）。層状ケイ酸塩 8 つは、違う大きさの超格子と 300 K の MD のあとの構造で全原子の型と電荷が正解と一致。参照を 1 つ抜くと金雲母だけ 16 原子（obts → ob）を間違える。
- 白雲母（ClayFF）：面内の弾性定数は実験より 30–50 % 大きく、層に垂直な引張は実験に近く、層をまたぐせん断は小さい。引張強度は x 14.7、y 13.4、z 3.5 GPa（破断ひずみ 0.102、0.079、0.103）で、ClayFF による先行研究（Teich-McGoldrick ら 2012、DelRio ら 2026）と一致。
- 拡散係数：水の O と Na⁺ は Greathouse ら（2016）の 0.75–0.99 倍、活性化エネルギーは誤差の範囲で一致。
