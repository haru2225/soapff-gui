# 例：ライブラリの 13 の構造の NVT と NPT

同梱の SOAP ライブラリ（`src/soapff_gui/data/library_clayff_bulk.npz`）を作った 13 の構造それぞれについて、構造ファイルから soapff-gui で力場を作り、NVT と NPT を回したものです。

| フォルダ | 鉱物 |
|---|---|
| `brucite/` | ブルーサイト（Mg(OH)₂） |
| `kaolinite/` | カオリナイト |
| `lizardite/` | リザーダイト |
| `montmorillonite_2W/` | Na モンモリロナイト 2 水和層（JAEA の構造、水と Na⁺ を含む） |
| `montmorillonite_dry_Na/` | 乾いた Na モンモリロナイト（Cygan ら 2004 の組成 Na₃(Si₃₁Al)(Al₁₄Mg₂)O₈₀(OH)₁₆、水なし。`examples/montmorillonite_swelling/` で作ったもの） |
| `muscovite_1M/`、`muscovite_2M1/` | 白雲母 1M、2M1 |
| `nanpingite/` | nanpingite（Cs 雲母） |
| `phlogopite_1M/` | 金雲母 1M |
| `pyrophyllite/`、`pyrophyllite_2/` | パイロフィライト（2 つの出典） |
| `quartz/` | 石英 |
| `talc/` | タルク |

各フォルダの `<名前>.extxyz` が構造ファイル（元素・座標・箱だけ。ライブラリの参照の LAMMPS データから、電荷と型を除いて書き出したもの）、`run/` が GUI の出力です。

## 手順（GUI）

1. 「新規計算」で `<名前>.extxyz` を選ぶ（同梱のライブラリ。力場と最小化までの `in.lammps` ができる）
2. 一覧から次を加える（温度 300 K、圧力 1 atm、時間刻み 0.5 fs）
   - ＋ その他 › soapff の緩和（初速度 → 緩和 NVT 10 ps）
   - ＋ その他 › トラジェクトリ（1 ps ごと）
   - ＋ fix › NVT、`run 100000`（50 ps）、`unfix`
   - ＋ fix › NPT（箱の傾きも含めて全部：結晶）または NPT（面に垂直な方向だけ：montmorillonite_2W）、`run 100000`（50 ps）、`unfix`
   - ＋ その他 › 最後の構造を書き出す
3. 「計算を開始」（このコンピュータ、4 並列（240 原子以下は 2 並列）、4–35 分）。終わると自動で解析（`run/analysis/`）

`run/` には、入力（`in.lammps`、`system.data`、`system.settings`、`charges.csv`、`report.json`）、出力（`log.lammps`、`final.data`）、解析（`analysis/`）があります。軌道（`traj.lammpstrj`）は大きいので入れていません。

## 結果

`python summarize.py` で `summary.csv` と下の表が出ます。温度・圧力・体積・密度は各 50 ps の後半の平均、体積の変化は NPT の後半の平均の、入力の箱に対する変化、PE の変化は NPT の後半の前と後ろの 12.5 ps の平均の差です。

| 鉱物 | 原子数 | fallback | 型の食い違い | T NVT (K) | T NPT (K) | P NPT (atm) | 体積の変化 (%) | 密度 (g/cm³) | PE の変化 (kcal/mol/原子) |
|---|---|---|---|---|---|---|---|---|---|
| brucite | 240 | 0 | 0 | 299.8 | 304.1 | 84 | +9.17 | 2.177 | +0.0051 |
| kaolinite | 408 | 0 | 0 | 300.9 | 299.1 | 263 | +0.67 | 2.591 | -0.0197 |
| lizardite | 324 | 0 | 0 | 301.3 | 297.9 | 25 | +3.58 | 2.492 | +0.0128 |
| montmorillonite_2W | 2502 | 0 | 0 | 299.6 | 301.6 | 5825 | +19.88 | 2.132 | +0.0006 |
| montmorillonite_dry_Na | 163 | 0 | 0 | 300.4 | 300.3 | 8591 | +0.61 | 2.815 | -0.0116 |
| muscovite_1M | 168 | 0 | 0 | 297.1 | 298.8 | -513 | -7.09 | 2.916 | -0.0048 |
| muscovite_2M1 | 336 | 0 | 0 | 300.4 | 301.0 | 310 | -1.96 | 2.904 | -0.0252 |
| nanpingite | 336 | 0 | 0 | 299.6 | 298.7 | -35 | -4.49 | 3.417 | -0.0031 |
| phlogopite_1M | 176 | 0 | 0 | 302.0 | 301.6 | 609 | +0.71 | 2.769 | +0.0112 |
| pyrophyllite | 320 | 0 | 0 | 301.3 | 297.6 | 353 | +1.96 | 2.761 | -0.0001 |
| pyrophyllite_2 | 320 | 0 | 0 | 303.3 | 301.2 | -522 | +2.09 | 2.757 | -0.0057 |
| quartz | 162 | 0 | 0 | 299.6 | 298.2 | 1381 | +4.74 | 2.530 | -0.0133 |
| talc | 378 | 0 | 0 | 299.5 | 303.8 | 91 | +2.85 | 2.699 | -0.0039 |

- すべての鉱物で、全原子が SOAP で決まり（fallback 0）、SOAP の型と結合のつながりの型の食い違いは 0 でした。
- どの鉱物も、温度は 300 K の ±5 K に収まり、NPT の後半でポテンシャルエネルギーはほぼ一定でした（0.03 kcal/mol/原子 以内）。
- NPT の体積の変化は、入力の構造（実験の格子や DFT の構造）と ClayFF の平衡の格子の違いです。ブルーサイト（+9 %）、白雲母 1M（−7 %）は大きく変わります（論文の 300 K の MD と同じ傾向）。
- montmorillonite_2W は、入力の構造の層間が狭く水が粘土に重なっているため、NPT で層間が広がり体積が 20 % 増えます。montmorillonite_2W と montmorillonite_dry_Na は x、y の箱を固定（面に垂直な方向だけの NPT）しているので、圧力は 1 atm になりません（乾いたものは 0 K で緩和した箱のまま）。
- アルバイトは ClayFF の対象外なので、ライブラリと例から外しました。
- 圧力の平均（P NPT）は、50 ps では ±1500 atm ほどゆらぎます（小さい結晶ほど大きい）。
