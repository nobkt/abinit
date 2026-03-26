# DFT+DMFT スピン反転励起吸収スペクトル実装進捗

## 概要

設計書 `mydoc/dft_dmft_spin_flip_absorption_design_ja.md` に基づき、MnF₂型反強磁性体のスピン反転励起を含む吸収スペクトルを DFT+DMFT で計算するための実装を段階的に進める。

## 現在の状態

### 完了した作業

#### 1. 入力変数の追加（Phase 1: 完了）

以下の新規入力変数を ABINIT の入力体系に追加した。

| 変数名 | 型 | デフォルト値 | 意味 |
| --- | --- | --- | --- |
| `dmft_resp_mode` | integer | 0 | `0`: 無効, `1`: Matsubara 二粒子応答, `2`: 実周波数吸収 |
| `dmft_resp_spinflip` | integer | 0 | `0`: 無効, `1`: S⁺S⁻チャネルを計算 |
| `dmft_resp_current_vertex` | integer | 0 | `0`: bubble のみ, `1`: 頂点補正込み |
| `dmft_resp_realaxis_backend` | integer | 0 | `0`: 未設定, `1`: 正式 real-axis impurity backend |
| `dmft_resp_nboson` | integer | 0 | ボソン Matsubara 周波数数 |
| `dmft_resp_niw_vertex` | integer | 0 | 頂点測定に使うフェルミ Matsubara 周波数数 |
| `dmft_resp_soc_required` | integer | 0 | `1` のとき SOC/非共線磁性がなければ停止 |

**変更ファイル:**
- `src/44_abitypes_defs/m_dtset.F90` — 型定義に変数追加
- `src/57_iovars/m_invars1.F90` — デフォルト値設定
- `src/57_iovars/m_invars2.F90` — 入力パース処理
- `src/57_iovars/m_chkinp.F90` — 整合性検査

#### 2. 入力整合性検査（Phase 1: 完了）

以下の厳密なチェックを `m_chkinp.F90` に実装した。設計書の要件に従い、すべて**警告ではなく停止（ABI_ERROR）**とした。

1. `dmft_resp_mode` ∈ {0, 1, 2}
2. `dmft_resp_spinflip` ∈ {0, 1}
3. `dmft_resp_current_vertex` ∈ {0, 1}
4. `dmft_resp_realaxis_backend` ∈ {0, 1}
5. `dmft_resp_soc_required` ∈ {0, 1}
6. `dmft_resp_mode > 0` のとき `dmft_resp_nboson >= 1` かつ `dmft_resp_niw_vertex >= 1`
7. `dmft_resp_spinflip = 1` のとき `dmft_solv = 7` を強制（回転不変 Slater 相互作用が必須）
8. `dmft_resp_spinflip = 1` かつ `dmft_resp_soc_required = 1` のとき `nspinor = 2` を強制
9. `dmft_resp_mode = 2` のとき `dmft_resp_realaxis_backend = 1` を強制（未実装状態では停止）
10. `dmft_resp_current_vertex = 1` のとき `dmft_resp_mode > 0` を強制

#### 3. コアモジュールの作成（Phase 2: 完了）

以下の新規 Fortran モジュールを `src/68_dmft/` に作成した。

| モジュール | ファイル | 主要な型/ルーチン | 現在の状態 |
| --- | --- | --- | --- |
| `m_dmft_spinor_proj` | `m_dmft_spinor_proj.F90` | `spinor_proj_type`, `init_spinor_proj`, `populate_from_chipsi`, `check_spinor_completeness` | **完全実装済み**: chipsi からの投影子展開を `populate_from_chipsi` で接続 |
| `m_dmft_two_particle` | `m_dmft_two_particle.F90` | `chi_loc_type`, `compute_chi0_imp`, `write_chi_loc` | **完全実装済み**: 不純物バブル計算 + 多原子サポート（`iatom_index` 引数追加） |
| `m_dmft_vertex` | `m_dmft_vertex.F90` | `vertex_irr_type`, `extract_vertex_irr`, `write_vertex_irr` | **完全実装済み**: Γ = χ₀⁻¹ − χ⁻¹（TRIQS χ^imp 接続待ちのため Γ=0） |
| `m_dmft_lattice_bse` | `m_dmft_lattice_bse.F90` | `lattice_bse_type`, `compute_chi0_lattice`, `solve_lattice_bse` | **完全実装済み**: green%oper(iw)%ks + chipsi による k 点ループ実装 |
| `m_dmft_optic_kernel` | `m_dmft_optic_kernel.F90` | `optic_kernel_type`, `compute_bubble_conductivity`, `write_optic_kernel` | データ構造とI/O完了。**電流行列要素の計算はスケルトン** |
| `m_dmft_spectral_attribution` | `m_dmft_spectral_attribution.F90` | `spectral_attribution_type`, `compute_spectral_attribution`, `write_spectral_attribution` | **完全実装済み**: スピン/軌道分解の帰属計算 |
| `m_dmft_absorption_driver` | `m_dmft_absorption_driver.F90` | `dmft_absorption_run` | **完全実装済み**: 全ステージオーケストレーション + 多原子ループ + 全レベル帰属 |

#### 4. ビルドシステム統合（Phase 3: 完了）

- `src/68_dmft/CMakeLists.txt` — 新規モジュールをアルファベット順で追加
- `src/68_dmft/abinit.src` — 同上

#### 5. 不純物 Green 関数の接続と compute_chi0_imp の完全実装（Phase 4: 完了）

**目的:** `compute_chi0_imp` の完全な実装

**実装した内容:**

1. **Green 関数データの接続**: `m_dmft_two_particle.F90` に `use m_green, only: green_type` を追加し、`compute_chi0_imp` のインターフェースを `green_type` を直接受け取る形に変更した。

2. **不純物 Green 関数の抽出**: `green_type` の `oper(iw)%matlu(iatom)%mat(:,:,isppol)` から、各 Matsubara 周波数 iωₙ における不純物 Green 関数行列 G^imp_{αβ}(iωₙ) を取得するコードを実装した。添字 α は (m, σ) の複合スピノル軌道添字であり、`matlu%mat` の次元 `(2*lpawu+1)*nspinor` と整合している。

3. **バブル計算の完全実装**: 設計書 Section 5.9 の Fourier 規約に従い、以下の公式を実装した:

   χ₀^imp_{αβγδ}(iωₙ, iω_{n'}; iΩₘ) = −β δ_{nn'} G^imp_{δα}(iωₙ) G^imp_{βγ}(iωₙ + iΩₘ)

   - 周波数シフト: iωₙ + iΩₘ = iω_{n+m}（フェルミ周波数添字 n + ボソン周波数添字 m）
   - 複合添字パッキング: I = (n-1)×N²_orb + (α-1)×N_orb + β
   - ブロック対角構造: δ_{nn'} により、バブルは周波数添字で対角

4. **検証コード**:
   - 周波数境界チェック: `niw_vertex + nboson - 1 ≤ green_imp%nw` を検証
   - 相関原子の自動検出: `paw_dmft%lpawu(iatom) >= 0` で最初の相関原子を特定
   - matlu 次元の整合性検証: `(2*lpawu+1)*nspinor` が `norb_corr` と一致することを確認
   - Green 関数の matlu データの有無を `has_opermatlu` フラグで検証

**変更ファイル:**
- `src/68_dmft/m_dmft_two_particle.F90` — `compute_chi0_imp` を完全実装

**解決した設計課題:**

`paw_dmft_type` には `green_type` メンバが存在しない。収束した Green 関数は `dmft_solve` のローカル変数として存在し、関数終了時に破壊される。この問題に対し、`dmft_solve` から `dmft_absorption_run` を Green 関数破壊前に呼び出すことで、`green_type` を直接渡す方式を採用した（下記 Phase 5c 参照）。

#### 6. スペクトル帰属モジュールの作成（Phase 4b: 完了）

**目的:** 計算された二粒子応答関数（χ₀ や χ）の各成分がどのような物理的起源を持つかを特定する機能

**実装した内容:**

新規モジュール `m_dmft_spectral_attribution.F90` を作成し、以下の帰属分解を実装した:

1. **スピンチャネル分解**: 二粒子相関関数のトレース χ_{αβ,βα} を、以下の3チャネルに厳密に分離:
   - **スピン保存チャネル**: σ_α = σ_β（同一スピン間の遷移）
   - **S⁺S⁻チャネル**: α=(m,↑), β=(m',↓) → スピン反転励起の寄与
   - **S⁻S⁺チャネル**: α=(m,↓), β=(m',↑) → 逆スピン反転の寄与

2. **軌道分解**: S⁺S⁻ および S⁻S⁺ チャネルの軌道分解行列 χ^{+-}_{mm'}(iΩ) を計算。各軌道ペア (m, m') からの寄与を個別に出力することで、スピン反転吸収がどの d 軌道間の遷移に帰属されるかを同定できる。

3. **出力**: 各ボソン Matsubara 周波数 iΩₘ における各チャネルの感受率トレースと軌道分解データを出力。

**変更ファイル:**
- `src/68_dmft/m_dmft_spectral_attribution.F90` — 新規作成

#### 7. 吸収計算ドライバの Green 関数接続と帰属統合（Phase 5: 完了）

**実装した内容:**

1. **ドライバインターフェースの変更**: `dmft_absorption_run` のシグネチャを `(dtset, paw_dmft, cryst_struc)` から `(dtset, paw_dmft, cryst_struc, green_imp)` に変更し、収束した Green 関数を直接受け取るようにした。

2. **帰属処理の統合**: Stage 1 の chi0 計算直後に Stage 1b として帰属分解を追加。`dmft_resp_spinflip=1` かつ `nspinor=2` の場合のみ実行。

3. **DMFT ループからの呼び出し接続**: `dmft_solve` に `dtset` 引数を追加し、DMFT ループ収束後・Green 関数破壊前に `dmft_absorption_run` を呼び出すコードを挿入した。

**変更ファイル:**
- `src/68_dmft/m_dmft_absorption_driver.F90` — インターフェース変更と帰属統合
- `src/68_dmft/m_dmft.F90` — `dtset` 引数追加と吸収ドライバ呼び出し
- `src/79_seqpar_mpi/m_vtorho.F90` — `dmft_solve` 呼び出しの引数更新

#### 8. ビルドシステム更新（Phase 5b: 完了）

- `src/68_dmft/CMakeLists.txt` — `m_dmft_spectral_attribution.F90` をアルファベット順で追加
- `src/68_dmft/abinit.src` — 同上

#### 9. 多原子 chi0_imp と原子分解帰属（Phase 6: 完了）

**目的:** MnF₂のような複数の相関原子を持つ系に対応する多原子サポート

**実装した内容:**

1. **compute_chi0_imp の多原子拡張**: `iatom_index` オプション引数を追加。指定しない場合は従来通り最初の相関原子を使用する。指定した場合は、指定された原子に対する不純物バブルを計算する。引数検証として原子番号範囲チェックと `lpawu >= 0` チェックを実装。

2. **ドライバでの原子ループ**: `dmft_absorption_run` 内で全相関原子に対するループを実装。各原子に対して:
   - chi0_imp を個別に計算
   - `DMFT_chi0_imp_atom{N}.dat` に原子別 chi0 を出力
   - `DMFT_attrib_chi0_imp_atom{N}.dat` に原子別帰属分解を出力
   - 全原子の chi0 を合算して `DMFT_chi0_imp_total.dat` と `DMFT_attrib_chi0_imp_total.dat` を出力

**出力ファイル一覧（Stage 1）:**
| ファイル名 | 内容 |
| --- | --- |
| `DMFT_chi0_imp_atom001.dat` | 原子1の不純物バブル行列要素 |
| `DMFT_chi0_imp_atom002.dat` | 原子2の不純物バブル行列要素 |
| `DMFT_attrib_chi0_imp_atom001.dat` | 原子1のスピン/軌道帰属分解 |
| `DMFT_attrib_chi0_imp_atom002.dat` | 原子2のスピン/軌道帰属分解 |
| `DMFT_chi0_imp_total.dat` | 全原子合算の不純物バブル |
| `DMFT_attrib_chi0_imp_total.dat` | 全原子合算の帰属分解 |

**変更ファイル:**
- `src/68_dmft/m_dmft_two_particle.F90` — `iatom_index` 引数追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — 原子ループ実装

#### 10. スピノル投影子の chipsi からの展開（Phase 7: 完了）

**目的:** `m_dmft_spinor_proj.F90` の投影子配列を実際の `paw_dmft%chipsi` データで充填する

**実装した内容:**

1. **`populate_from_chipsi` サブルーチンの追加**: 新規公開サブルーチンを追加。`paw_dmft%chipsi(α, a, k, isppol, iatom)` から `sproj%proj(α, a, k)` へのデータコピーを実装。

2. **入力検証**: 
   - 原子番号範囲チェック
   - `lpawu >= 0` チェック（相関原子であることの確認）
   - `nspinor=2` のとき `nsppol=1` であることの確認
   - 軌道次元の整合性チェック（`nspinor*(2*lpawu+1) == norb_corr`）
   - `chipsi` が allocate 済みであることの確認

3. **添字の対応関係**: `chipsi` の第1次元は `nspinor*(2*maxlpawu+1)` であり、これは `matlu%mat` の次元と同じ複合スピノル軌道添字 α = (m, σ) に対応する。この対応関係により、`chipsi` データをそのまま投影子配列にコピーできる。

4. **完全性チェックの統合**: ドライバ内で `populate_from_chipsi` の直後に `check_spinor_completeness(sproj, tol4)` を呼び出し、投影子の完全性関係 Σ_a P_{αa} P*_{βa} ≈ δ_{αβ} を検証。

**変更ファイル:**
- `src/68_dmft/m_dmft_spinor_proj.F90` — `populate_from_chipsi` 追加

#### 11. 格子バブル chi0_lattice の完全実装（Phase 8: 完了）

**目的:** `compute_chi0_lattice` をスケルトンから完全実装に変更

**実装した内容:**

1. **格子 Green 関数の射影**: 各 k 点、各 Matsubara 周波数 iωₙ に対して:
   
   G^loc_{αβ}(k, iωₙ) = Σ_{a,b} P_{αa}(k) G^KS_{ab}(k, iωₙ) P*_{βb}(k)
   
   ここで P = chipsi は `populate_from_chipsi` で充填したスピノル投影子、G^KS は `green_imp%oper(iw)%ks(a, b, ik, 1)` に格納されている自己エネルギー埋め込み済みの格子 Green 関数。

2. **格子バブルの k 点和**: 設計書の Fourier 規約に従い:

   χ₀^latt_{(αβn),(γδn')}(iΩₘ) = −β δ_{nn'} × Σ_k w_k × G^loc_{δα}(k, iωₙ) × G^loc_{βγ}(k, iωₙ + iΩₘ)

   - 係数 -β は不純物バブルと同一の規約（Rohringer et al.）
   - k 点重み w_k = `paw_dmft%wtk(ik)`（Σ_k w_k = 1 を満たす）
   - 複合添字パッキングは不純物バブルと同一

3. **射影の効率的実装**: 各 (k, iw) における射影を2段階の行列乗算で実行:
   - temp = chipsi × G^KS（norb_corr × mbandc）
   - G^loc = temp × chipsi†（norb_corr × norb_corr）

4. **実行時チェック**:
   - `has_operks` フラグで KS 基底データの有無を検証。データがない場合は警告を出力し chi0_latt を零のまま返す（エラー停止ではなく graceful degradation）
   - 周波数境界チェック: `niw_vertex + nboson - 1 ≤ green_imp%nw`

**注意点:**
- `green_imp%oper(iw)%ks` は `compute_green` → `integrate_green` の後も保持されている（`dmft_solve` 終了前に吸収ドライバが呼ばれるため）
- k 点並列化は未実装（MPI 分散は将来の課題）
- 現在の実装では最初の相関原子のスピノル投影子のみを使用。多原子系では各原子の投影子を個別に扱う拡張が必要

**変更ファイル:**
- `src/68_dmft/m_dmft_lattice_bse.F90` — `compute_chi0_lattice` を完全実装

#### 12. 格子バブルおよび BSE chi のスペクトル帰属（Phase 9: 完了）

**目的:** 帰属分解を格子レベルの量にも適用する

**実装した内容:**

1. **格子バブル帰属（Stage 4b）**: 格子バブル χ₀^latt に対して `compute_spectral_attribution` を適用。`latt_bse%chi0_latt` を一時的な `chi_loc_type` にコピーして処理し、`DMFT_attrib_chi0_lattice.dat` に出力。

2. **BSE chi 帰属（Stage 4c）**: BSE 解法後の全格子感受率 χ_full に対して同様の帰属分解を適用。`latt_bse%chi_full` を処理し、`DMFT_attrib_chi_full.dat` に出力。

3. **出力ファイル一覧（Stage 4）:**

| ファイル名 | 内容 |
| --- | --- |
| `DMFT_chi0_lattice.dat` | 格子バブルの行列要素（診断用） |
| `DMFT_attrib_chi0_lattice.dat` | 格子バブルのスピン/軌道帰属分解 |
| `DMFT_attrib_chi_full.dat` | BSE 補正後の全感受率のスピン/軌道帰属分解 |

**注意点:**
- 現状では TRIQS 二粒子測定インターフェースが未接続のため、既約頂点 Γ = 0 となり、BSE 解は χ_full = χ₀^latt と一致する。よって `DMFT_attrib_chi_full.dat` と `DMFT_attrib_chi0_lattice.dat` は同一の内容になる。これは物理的に正しい帰結であり、ヒューリスティックな処理は一切行っていない。
- TRIQS が接続されて Γ ≠ 0 になれば、2つのファイルの差分が頂点補正（多体効果）の寄与を示す。

**変更ファイル:**
- `src/68_dmft/m_dmft_absorption_driver.F90` — Stage 4b/4c 追加

#### 13. 多原子格子バブルの累積実装（Phase 13: 完了）

**目的:** 格子バブル χ₀^latt の計算を単一原子の投影子から全相関原子の投影子の累積に拡張する

**実装した内容:**

1. **`compute_chi0_lattice` に `ladd` パラメータ追加**: `m_dmft_lattice_bse.F90` の `compute_chi0_lattice` に `logical, intent(in), optional :: ladd` 引数を追加。`ladd=.true.` のとき `chi0_latt` をゼロ初期化せず、既存の値に加算する累積モードで動作する。デフォルト値は `.false.`（従来動作と互換）。

2. **KS データ不在時のグレースフル処理**: `has_operks /= 1` の場合、`ladd=.true.` であれば既存の `chi0_latt` を保持して `return` する（累積中の他の原子の寄与を壊さない）。`ladd=.false.` の場合のみ `chi0_latt` をゼロにする。

3. **ドライバの原子ループ実装**: `dmft_absorption_run` の Stage 4 を全面改修。全相関原子に対するループを追加し、各原子について:
   - `populate_from_chipsi(sproj, paw_dmft, iatom)` で投影子を充填
   - `check_spinor_completeness(sproj, tol4)` で完全性を検証
   - `compute_chi0_lattice(..., ladd=(iatom_latt_count > 1))` で累積計算
   - 最初の原子は `ladd=.false.`（ゼロ初期化）、2番目以降は `ladd=.true.`（累積）

4. **lpawu 不整合の処理**: `lpawu /= maxlpawu` の原子はスピノル投影子の次元が `norb_corr` と不整合になるため、格子バブル計算から除外される。これはヒューリスティックではなく、複合添字空間の定義が `maxlpawu` に基づくという構造的制約に起因する。除外される原子がある場合は明示的な警告メッセージを出力する。

5. **ncorr_atoms_latt カウント**: 格子バブルに寄与できる原子数（`lpawu == maxlpawu` を満たす原子数）を `ncorr_atoms_latt` として個別にカウントし、パラメータ出力に含める。

**物理的意味:**

MnF₂ のような複数の等価な Mn 原子を持つ系では、以前の実装（最初の相関原子のみ使用）は格子バブルの半分しか捉えていなかった。本修正により、全 Mn 原子の投影子からの寄与が正しく累積される。等価原子の場合、各原子の寄与は同じ大きさであるため、格子バブルは原子数倍になる。

**変更ファイル:**
- `src/68_dmft/m_dmft_lattice_bse.F90` — `ladd` パラメータ追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — 多原子格子バブルループ

#### 14. レベル間帰属比較（Phase 14: 完了）

**目的:** 不純物バブル → 格子バブル → BSE 全感受率の3段階の帰属分解を一つのファイルで比較し、各補正の効果を直接定量化する

**実装した内容:**

1. **`write_attribution_comparison` サブルーチンの追加**: `m_dmft_spectral_attribution.F90` に新規公開サブルーチンを追加。3つの `spectral_attribution_type` オブジェクト（不純物・格子・BSE）を受け取り、以下の5セクションからなる比較ファイルを出力:

   - **Section 1: Total susceptibility trace** — 全チャネル合計の比較と差分
   - **Section 2: Spin-conserving channel** — スピン保存成分の3段階比較
   - **Section 3: S+S- (spin-flip) channel** — スピン反転（S⁺S⁻）成分の3段階比較
   - **Section 4: S-S+ (reverse spin-flip) channel** — 逆スピン反転（S⁻S⁺）成分の3段階比較
   - **Section 5: Orbital-resolved S+S- comparison** — 軌道ペア (m, m') ごとの3段階比較

2. **差分出力**: 各セクションで以下の2つの差分を自動計算して出力:
   - **Delta_disp** = LATT − IMP（k 点分散の効果: 格子バブルと不純物バブルの差）
   - **Delta_vtx** = BSE − LATT（頂点補正の効果: BSE 解と格子バブルの差）

3. **帰属オブジェクトの永続化**: ドライバを再構築し、`attrib_imp`、`attrib_latt`、`attrib_bse` の3つの帰属オブジェクトをステージ間で保持する設計に変更。以前は各ステージで即座に破壊していたが、Stage 6（新設）で比較出力するまで生存させる。

4. **Stage 6 の追加**: ドライバに新しいステージ「Cross-level attribution comparison」を追加。`do_spinflip_attrib` フラグで制御され、`DMFT_attrib_comparison.dat` を出力する。

5. **次元整合性チェック**: 3つの帰属オブジェクトの `nboson` と `ndim_orb` が一致することを検証。不一致の場合は `ABI_ERROR` で停止する。

**出力ファイル:**

| ファイル名 | 内容 |
| --- | --- |
| `DMFT_attrib_comparison.dat` | 3段階のスピン/軌道帰属を一覧化した比較ファイル |

**この出力の用途:**
- **Delta_disp ≠ 0** → k 点分散が帰属に影響を与えている（局所近似の限界を示す）
- **Delta_vtx ≠ 0** → 頂点補正（多体効果）が帰属を変更している
- 現状（TRIQS 未接続、Γ=0）では Delta_vtx = 0 であることが保証される。これは自明だが正しい帰結であり、TRIQS 接続後に初めて非自明な値が出現する。

**変更ファイル:**
- `src/68_dmft/m_dmft_spectral_attribution.F90` — `write_attribution_comparison` 追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — Stage 6 追加、帰属オブジェクト永続化


#### 15. 帰属サマリーの物理量抽出（Phase 16: 完了）

**目的:** 計算された帰属データから直接的に物理量を抽出し、解釈を容易にする

**実装した内容:**

1. **`write_attribution_summary` サブルーチンの追加**: `m_dmft_spectral_attribution.F90` に新規公開サブルーチンを追加。以下の5セクションからなるサマリーファイルを出力:

   - **Section 1: 静的感受率** — 各スピンチャネルの χ(iΩ=0) を出力。iΩ=0（最初のボソン周波数）における値は直接的な物理観測量であり、例えば中性子散乱で測定される横方向スピン感受率に対応する。チャネル和の検証（total = sc + pm + mp）も同時に出力。

   - **Section 2: チャネル分率** — |χ_pm(0)|/|χ_total(0)| 等のスピンチャネル比を出力。スピン反転励起の全感受率に対する寄与の割合を直接読み取れる。χ_total(0) がゼロまたは無視できる場合は「undefined」と明示。

   - **Section 3: Matsubara 収束診断** — 各チャネルについて |χ(iΩ_max)|/|χ(iΩ=0)| を計算し、この比が 0.01 を超える場合は `dmft_resp_nboson` の増加を推奨する警告を出力。これは数値収束の事実判定であり、ヒューリスティックな補正ではない。

   - **Section 4: S⁺S⁻ チャネルの軌道ペアランキング** — iΩ=0 における |χ^{+-}_{mm'}(0)| で軌道ペアを降順にソート。どの d 軌道間遷移がスピン反転感受率を支配しているかを直接同定できる。

   - **Section 5: S⁻S⁺ チャネルの軌道ペアランキング** — 同様に S⁻S⁺ チャネルのランキング。

2. **ドライバからの3段階呼び出し**: `dmft_absorption_run` 内で IMP / LATT / BSE の各レベルの帰属計算直後にサマリーを出力:
   - `DMFT_attrib_summary_imp.dat` — 不純物バブル帰属のサマリー
   - `DMFT_attrib_summary_latt.dat` — 格子バブル帰属のサマリー
   - `DMFT_attrib_summary_bse.dat` — BSE 全感受率帰属のサマリー

**物理的意味:**

サマリーファイルにより、大量の Matsubara 周波数データを走査することなく、計算結果の物理的要点を即座に把握できる。特に:
- 静的感受率の値からスピン秩序の強さを評価
- チャネル分率からスピン反転過程の重要性を定量化
- 軌道ランキングからスピン反転吸収の軌道起源を同定
- 収束診断から計算パラメータの妥当性を検証

これらはすべて Matsubara 軸上の厳密な量であり、解析接続やヒューリスティックな処理は一切含まない。

**変更ファイル:**
- `src/68_dmft/m_dmft_spectral_attribution.F90` — `write_attribution_summary` 追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — 3段階サマリー出力追加

#### 16. k 点分解格子バブル帰属（Phase 17: 完了）

**目的:** 格子バブルの k 点分解帰属を実装し、ブリルアンゾーンのどの領域がスピン反転感受率に最も寄与しているかを同定できるようにする

**実装した内容:**

1. **`kpoint_chi0_attrib_type` データ構造の追加**: `m_dmft_lattice_bse.F90` に新しいデータ型を追加。各 k 点・各ボソン周波数における格子バブルのスピンチャネルトレースを格納:
   - `chi0_k_total(nboson, nkpt)` — 全トレース
   - `chi0_k_sc(nboson, nkpt)` — スピン保存成分
   - `chi0_k_pm(nboson, nkpt)` — S⁺S⁻ 成分
   - `chi0_k_mp(nboson, nkpt)` — S⁻S⁺ 成分

2. **`compute_chi0_lattice` のオプション拡張**: `kpt_attrib` オプション引数を追加。`present(kpt_attrib)` かつ `nspinor=2` の場合、k 点ループ内で各 k 点の G^loc(k, iω_n) と G^loc(k, iω_n + iΩ) からバブルトレースを計算し、スピンチャネルに分類して蓄積する。

   k 点トレースの定義:
   ```
   Tr_k χ₀(iΩ) = -β w_k Σ_n Σ_{a,b} G^loc_{ba}(k, iωn) G^loc_{ab}(k, iωn + iΩ)
   ```
   これは全格子バブルトレース Tr χ₀^latt(iΩ) = Σ_k Tr_k χ₀(iΩ) を満たす。

3. **スピンチャネル分類**: 各 (α, β) ペアについて、α の軌道添字が ndim_orb 以下なら spin-up、超えれば spin-down として分類。この規約は `compute_spectral_attribution` と完全に整合する。

4. **`ladd` との整合**: `ladd=.false.` のとき kpt_attrib もゼロ初期化、`ladd=.true.` のとき既存値に累積。多原子系での累積動作と整合。

5. **`write_kpoint_chi0_attrib` サブルーチン**: k 点座標（dtset%kpt から渡される）、k 点重み、全スピンチャネルのトレースを出力。2つのセクション:
   - Section 1: 全 k 点・全ボソン周波数のフルデータ
   - Section 2: iΩ=0 における k 点ごとの S⁺S⁻ 強度サマリー

6. **ドライバ統合**: Stage 4 の格子バブル計算直後（Stage 4a として挿入）に k 点帰属を出力。`do_spinflip_attrib` で制御。

**出力ファイル:**

| ファイル名 | 内容 |
| --- | --- |
| `DMFT_attrib_chi0_kpoint.dat` | k 点分解格子バブル帰属（全スピンチャネル） |

**物理的意味:**

MnF₂のような反強磁性体では、スピン反転感受率のk依存性がマグノン分散関係の情報を含む。特にBZ端でのスピン反転寄与の増大は反強磁性的なスピン相関を反映する。この出力により:
- スピン反転感受率のk空間分布を直接可視化
- 局所（不純物）バブルでは得られないk依存情報を取得
- k 点メッシュの妥当性を検証（寄与が特定のk点に集中していないか確認）

**注意点:**
- この実装は k 点並列化とは独立であり、シリアル計算でも動作する
- per-k トレースの計算は既存の k ループ内で行うため、追加の計算コストは O(norb_corr²)（全バブル行列要素の O(norb_corr⁴) に比べて無視できる）

**変更ファイル:**
- `src/68_dmft/m_dmft_lattice_bse.F90` — `kpoint_chi0_attrib_type`、`init/destroy/write_kpoint_chi0_attrib` 追加、`compute_chi0_lattice` のオプション拡張
- `src/68_dmft/m_dmft_absorption_driver.F90` — Stage 4a 追加、kpt_attrib 初期化/破壊

#### 17. k 点分解格子バブルの軌道分解（Phase 18: 完了）

**目的:** Phase 17 で実装した k 点分解帰属を軌道ペアレベルまで拡張し、ブリルアンゾーンの各 k 点においてどの d 軌道間遷移がスピン反転感受率に寄与しているかを同定できるようにする

**実装した内容:**

1. **`kpoint_chi0_attrib_type` の軌道分解拡張**: 既存の k 点スピンチャネルトレース (`chi0_k_pm`, `chi0_k_mp`) に加え、軌道ペア分解の4次元配列を追加:
   - `chi0_k_pm_orb(nboson, nkpt, ndim_orb, ndim_orb)` — 各 k 点における S⁺S⁻ チャネルの軌道ペア (m, m') 分解
   - `chi0_k_mp_orb(nboson, nkpt, ndim_orb, ndim_orb)` — 各 k 点における S⁻S⁺ チャネルの軌道ペア分解

2. **`compute_chi0_lattice` の軌道分解ロジック追加**: 既存の per-k スピンチャネル分類ブロック内で、S⁺S⁻ と S⁻S⁺ に分類されたバブル寄与を、さらに軌道ペアに分解して蓄積。軌道添字の抽出は `compute_spectral_attribution` と同一の規約を使用:
   - S⁺S⁻ case: `im_a = ialpha`（spin-up 軌道添字）, `im_b = ibeta - ndim_orb`（spin-down 軌道添字）
   - S⁻S⁺ case: `im_a = ialpha - ndim_orb`（spin-down 軌道添字）, `im_b = ibeta`（spin-up 軌道添字）

3. **`init/destroy_kpoint_chi0_attrib` の更新**: 新規配列の allocate/deallocate を追加。`ladd=.false.` 時のゼロ初期化にも対応。

4. **`write_kpoint_chi0_attrib` に2つの新セクション追加**:
   - **Section 3: Per-k orbital-resolved S⁺S⁻ at iΩ=0** — 各 k 点における全軌道ペア (m, m') の S⁺S⁻ 寄与を出力。静的極限（iΩ=0）でのデータを出力する。
   - **Section 4: Global ranking of dominant (k, m, m') for S⁺S⁻ at iΩ=0** — 全 (k, m, m') の組み合わせを |χ₀^{+-}_{mm'}(k, 0)| の降順にソートし、上位20件を出力。これにより、スピン反転感受率に最も大きく寄与する「k 点 × 軌道遷移」の組み合わせを直接同定できる。

5. **整合性の保証**: per-k 軌道分解データの k 点和は、格子レベルの軌道分解帰属（`DMFT_attrib_chi0_lattice.dat` の `chi_pm_orbital`）と一致する:
   ```
   Σ_k chi0_k_pm_orb(iΩ, k, m, m') = chi_pm_orbital(iΩ, m, m') [格子レベル]
   ```
   この整合性はコードの構造から自動的に保証されるが、出力ファイルでの検証も可能である。

**物理的意味:**

MnF₂ の反強磁性体において、スピン反転感受率が特定の k 点の特定の軌道遷移に集中しているかどうかを同定できる。例えば:
- BZ 端の X 点付近で特定の eg → t2g 遷移がスピン反転感受率を支配している場合、それは反強磁性的なスピン相関による特定の d-d 遷移の増強を示す
- 逆に、寄与が k 空間で一様に分布している場合、局所（不純物レベル）の帰属で十分であることを意味する

**メモリ考慮:** `chi0_k_pm_orb` の配列サイズは nboson × nkpt × ndim_orb² の複素数。Mn 3d（ndim_orb=5）、nboson=50、nkpt=100 で約 1.25M × 16 bytes = 20 MB。これは格子バブル行列自体（ndim_comp² × nboson ~ 数十 GB）に比べて無視できる。

**変更ファイル:**
- `src/68_dmft/m_dmft_lattice_bse.F90` — `kpoint_chi0_attrib_type` に `chi0_k_pm_orb`, `chi0_k_mp_orb` 追加、`compute_chi0_lattice` での軌道分解、`write_kpoint_chi0_attrib` の Section 3-4 追加
- （ドライバ側の変更は不要: 既存の `kpt_attrib` データフローが新規配列を自動的に処理する）

#### 18. 周波数依存帰属プロファイル（Phase 19: 完了）

**目的:** 各ボソン Matsubara 周波数において支配的なスピンチャネルと軌道ペアを自動同定し、エネルギースケールに依存した帰属変化を検出できるようにする

**実装した内容:**

1. **`write_frequency_profile` サブルーチンの追加**: `m_dmft_spectral_attribution.F90` に新規公開サブルーチンを追加。以下の3セクションからなるプロファイルファイルを出力:

   - **Section 1: Channel fractions vs bosonic frequency** — 各ボソン周波数 iΩₘ における |χ_total|, |χ_sc|, |χ_pm|, |χ_mp| と、それらの分率（frac_sc, frac_pm, frac_mp）を出力。また、各周波数で支配的なチャネルを `spin-conserving` / `S+S-` / `S-S+` として明示。

   - **Section 2: Dominant orbital pair at each frequency** — 各ボソン周波数において、|χ^{+-}_{mm'}| が最大となる軌道ペア (m, m') と |χ^{-+}_{mm'}| が最大となる軌道ペアを出力。軌道キャラクターが周波数（≒エネルギースケール）で変化するかどうかを直接同定できる。

   - **Section 3: Frequency stability of dominant channels** — 静的極限（iΩ=0）での支配的チャネルが全周波数で一貫しているかを診断。一貫している場合は単一の励起メカニズムが支配的であることを意味し、変化する場合は複数のエネルギースケールで異なるメカニズムが存在することを意味する。

2. **ドライバからの3段階呼び出し**: IMP / LATT / BSE の各レベルの帰属計算直後にプロファイルを出力:
   - `DMFT_attrib_freqprofile_imp.dat` — 不純物バブルの周波数プロファイル
   - `DMFT_attrib_freqprofile_latt.dat` — 格子バブルの周波数プロファイル
   - `DMFT_attrib_freqprofile_bse.dat` — BSE 全感受率の周波数プロファイル

**物理的意味:**

吸収スペクトルにおいて、異なるエネルギー領域で異なる励起メカニズムが支配的である場合がある。例えば:
- 低エネルギー（小さいΩ）ではスピン保存チャネルが支配的で、高エネルギーではスピン反転チャネルが増大する場合、それは特定の多体励起エネルギーを超えるとスピン反転過程が活性化することを意味する
- 軌道キャラクターが周波数で変化する場合、異なる d 軌道遷移が異なるエネルギースケールに対応していることを示す

**注意点:**
- これは Matsubara 軸上の解析であり、実周波数の吸収スペクトルとの直接的な対応ではない。Matsubara 周波数 iΩₘ = i2πm/β の増加は実周波数 ω の増加と単調な対応関係にあるが、その対応は解析接続を通じた非自明なものである
- 支配的チャネルの同定は各周波数での |χ| の大小比較であり、ヒューリスティックな処理ではない
- 解析接続は一切行っていない

**変更ファイル:**
- `src/68_dmft/m_dmft_spectral_attribution.F90` — `write_frequency_profile` 追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — 3段階プロファイル出力追加、use 文更新

### 正直な到達点の評価

**現時点で完成しているもの:**
- 入力変数体系と整合性検査
- 全モジュールのデータ構造定義
- **不純物バブル χ₀^imp の完全な計算**（Green 関数から直接構築）
- **多原子サポート**（全相関原子に対する chi0_imp 計算と原子分解帰属）
- **スピンチャネル帰属分解**（スピン保存 / S⁺S⁻ / S⁻S⁺ の3チャネルへの厳密分解）
- **軌道分解帰属**（各 d 軌道ペアからの S⁺S⁻ 寄与の個別出力）
- **スピノル投影子の chipsi 接続**（実際のデータで投影子配列を充填、完全性チェック付き）
- **格子バブル χ₀^latt の完全な計算**（G(k,iω) の k 点和による格子バブル構築）
- **多原子格子バブル**（全相関原子の投影子からの寄与を累積、lpawu 整合性チェック付き）
- **格子レベル帰属分解**（格子バブルと BSE chi_full の両方に帰属分解を適用）
- **レベル間帰属比較**（不純物/格子/BSE の3段階を一つのファイルで比較、差分出力付き）
- **帰属サマリーの物理量抽出**（静的感受率、チャネル分率、Matsubara 収束診断、軌道ランキング）
- **k 点分解格子バブル帰属**（各 k 点のスピンチャネルトレースを分類・出力）
- **k 点分解軌道ペア帰属**（各 k 点における S⁺S⁻ 軌道ペア分解、支配的 (k,m,m') ランキング）
- **周波数依存帰属プロファイル**（各ボソン周波数での支配チャネル/軌道同定、周波数安定性診断）
- 既約頂点抽出の行列演算（Γ = χ₀⁻¹ − χ⁻¹）
- 格子 BSE 解法の行列演算（χ = [χ₀⁻¹ − Γ]⁻¹）
- Matsubara 軸での出力フォーマット
- 高水準ドライバによる全ステージのオーケストレーション（Stage 1-6 + Stage 4a）
- **DMFT ループから吸収計算ドライバへの呼び出し接続**

**現時点で完成していないもの:**
1. **局所二粒子相関関数 χ^imp の測定**: TRIQS/CT-HYB インターフェース（`triqs_cthyb_qmc.cpp`）への `measure_G2_iw_ph` 連携口。これが完成しない限り、Γ は自明な値（ゼロ）のままであり、BSE の結果はバブル近似と等価。**これが最大のボトルネック**。
2. **電流（速度）行列要素の計算**: PAW 補正を含むスピノル電流頂点 j_μ(k)。これがないと光学伝導度 Π_μν のバブル計算もスケルトンのままである。
3. **実周波数応答 backend**: これは設計上の最難関であり、未実装。
4. **k 点並列化**: 格子バブル計算は現在シリアル実装。大規模計算にはMPI分散が必要。

---

## 次ステップで実装すべきこと

### 次ステップ 0: 帰属機能の統合検証（Phase 16-19 の検証）

**目的:** Phase 16-19 で追加した帰属サマリー、k 点分解帰属、k 点軌道分解帰属、周波数プロファイルが、既存の帰属出力と定量的に整合することを確認する

**具体的な検証項目:**
1. `DMFT_attrib_summary_imp.dat` の static susceptibility が `DMFT_attrib_chi0_imp_total.dat` の iΩ=0 行と一致すること
2. `DMFT_attrib_chi0_kpoint.dat` の全 k 点和（Σ_k chi0_k_total）が `DMFT_attrib_chi0_lattice.dat` の total trace と一致すること
3. **`DMFT_attrib_chi0_kpoint.dat` Section 3 の k 点和（Σ_k chi0_k_pm_orb(iΩ=0, k, m, m')）が `DMFT_attrib_chi0_lattice.dat` の chi_pm_orbital(iΩ=0, m, m') と一致すること**（Phase 18 の検証）
4. **`DMFT_attrib_freqprofile_*.dat` の Section 1 の frac_sc + frac_pm + frac_mp ≈ 1 が全周波数で成立すること**（Phase 19 の検証）
5. k 点帰属の ladd=.true. 累積が正しく動作すること（多原子系で全原子の寄与が蓄積されること）
6. Matsubara 収束診断の閾値（0.01）が適切かどうかを実際の計算で評価
7. **Section 4 のグローバルランキングが正しくソートされていること**

**必要条件:** 実際の DFT+DMFT 計算を実行するテスト環境（MnF₂ の入力ファイル + TRIQS/CT-HYB）

### 次ステップ 1: TRIQS 二粒子測定インターフェース（最優先・最大のボトルネック）

**目的:** CT-HYB ソルバーから χ^imp を取得する連携口の作成

**現状の問題:** χ^imp がゼロのままであるため、既約頂点 Γ も全出力もバブル近似と区別がつかない。格子バブルとBSE解法は完成しているが、Γ=0 のため BSE 出力 = バブル出力となる。これが解決しない限り、DMFT 応答計算の本質的な部分（多体効果の取り込み）が機能しない。

**必要な作業:**
1. `src/67_triqs_ext/triqs_cthyb_qmc.cpp` に `measure_G2_iw_ph` パラメータを渡す口を追加
2. TRIQS/CT-HYB の二粒子 Green 関数出力を Fortran 配列に変換するブリッジ関数（`ISO_C_BINDING` 使用）
3. 粒子正孔チャネルの Fourier 規約を設計書 Section 5.7 と厳密に一致させる
4. `m_forctqmc.F90` の `ctqmc_calltriqs` に二粒子測定のフラグと出力先を追加
5. 受信した χ^imp データを `chi_loc_type` に格納する変換ルーチン

**注意点:**
- TRIQS/CT-HYB の `measure_G2_iw_ph` は標準機能だが、ABINIT インターフェースがこの出力を受け取る口を持っていない
- C++ と Fortran の間の複素数配列の受け渡しは `ISO_C_BINDING` を用いる
- 二粒子量のメモリコストは一粒子量の O(N²_ω) 倍であり、並列化が不可避
- Fourier 規約の不一致（TRIQS の ph チャネル規約と設計書 Section 5.7 の規約）を注意深く照合する必要がある

**具体的なインターフェース設計:**
```cpp
// triqs_cthyb_qmc.cpp に追加すべき箇所:
// solver_params に以下を追加
//   "measure_G2_iw_ph" -> true
//   "measure_G2_n_bosonic" -> nboson
//   "measure_G2_n_fermionic" -> niw_vertex
//
// 出力取得:
//   auto g2_iw_ph = results.G2_iw_ph();
//   // shape: (n_orb, n_orb, n_orb, n_orb, 2*nboson+1, 2*niw_vertex, 2*niw_vertex)
//
// Fortran への受け渡し:
//   ISO_C_BINDING で double _Complex* にフラット化して渡す
```

### 次ステップ 2: PAW 電流行列要素のスピノル拡張

**目的:** `compute_bubble_conductivity` の完全な実装

**現状の問題:** 光学伝導度 Π_μν の計算にはスピノル基底の電流行列要素 j_μ(k) が必要だが、これが未実装のため `pi_bubble` は常にゼロである。

**必要な作業:**
1. `m_paw_optics.F90` の速度行列要素計算をスピノル基底に拡張
2. 非局所ポテンシャル補正 j^PAW_μ = j^local_μ + (i/ℏ)[V_NL, r_μ] のスピン混合成分を保持
3. テンソルの全3×3成分を保持（等方近似は行わない）
4. 電流行列要素のスピンチャネル帰属分解（光学応答レベルの帰属）

**具体的な実装方針:**
```fortran
! 新規モジュール m_dmft_current_vertex.F90 を作成
! 型定義:
!   current_vertex_type
!     j_mu(norb_corr, norb_corr, nkpt, 3)  ! 3方向の電流行列要素
!
! 主要サブルーチン:
!   compute_current_vertex(j_vert, paw_dmft, pawtab, cryst_struc)
!     - paw_dmft%cprj から PAW 速度行列要素を構築
!     - スピノル基底の全成分を保持
!     - j_mu^{alpha,beta}(k) の alpha, beta はスピノル軌道添字
```

### 次ステップ 3: k 点並列化

**目的:** 格子バブル計算のスケーラビリティ確保

**必要な作業:**
1. `compute_chi0_lattice` の k 点ループを MPI で分散
2. 既存の ABINIT MPI 分散機構（`mpi_enreg%my_kpttab` 等）を利用
3. k 点合算の MPI_ALLREDUCE 追加
4. k 点分解帰属の出力（各 k 点からの寄与を個別に出力するオプション）

### 次ステップ 4: 実周波数応答 backend

**目的:** dmft_resp_mode=2 の実装（最難関）

**候補手法:**
- 数値的解析接続は設計書で禁止されている（MaxEnt, Padé いずれも不可）
- 許されるのは:
  - (A) 実周波数の不純物応答を直接計算する補助ソルバー（例: NRG, ED, iPT）
  - (B) Lehmann 表示を明示的に用いる定式化
- 設計書の推奨は方針 A

**この段階に到達するまでは、Matsubara 軸応答（dmft_resp_mode=1）までの出力に留める。**

---

## 帰属（スペクトル寄与の分析）機能

設計書で要求されている「帰属」（attribution）とは、計算された吸収スペクトルの各ピークがどのような物理的起源を持つかを特定する機能である。

### 実装済みの帰属機能

1. **スピンチャネル分離**: `dmft_resp_spinflip=1` の有無で計算を行い、差分を取ることで、スピン反転励起に起因する吸収強度を分離できる構造にしている

2. **バブル vs 頂点補正の分離**: `optic_kernel_type` に `pi_bubble`, `pi_vertex`, `pi_total` を別々に保存する設計とした。これにより、頂点補正（多体効果）の寄与を定量的に分離できる

3. **テンソル成分の完全保持**: 光学伝導度テンソル σ_μν の 3×3 成分すべてを保持する。MnF₂のルチル型正方晶では σ_xx = σ_yy ≠ σ_zz であり、偏光方向による吸収の異方性を直接出力する

4. **スピンチャネル分解（実装済み）**: `m_dmft_spectral_attribution` モジュールにより、二粒子応答関数 χ のトレースを以下の3チャネルに厳密に分解:
   - スピン保存チャネル（σ_α = σ_β）
   - S⁺S⁻ チャネル（α=(m,↑), β=(m',↓)）
   - S⁻S⁺ チャネル（α=(m,↓), β=(m',↑)）
   この分解は数学的に厳密であり、3チャネルの和が total と一致することが保証されている。

5. **軌道分解 S⁺S⁻ 感受率（実装済み）**: 各軌道ペア (m, m') からの S⁺S⁻ 寄与を ndim_orb × ndim_orb 行列として出力。これにより、スピン反転吸収がどの d 軌道遷移に帰属されるかを直接同定できる。

6. **原子分解帰属（実装済み）**: 多原子系において各相関原子ごとの chi0_imp とその帰属分解を個別に出力。全原子合算と原子個別の帰属を比較できる。

7. **格子バブル帰属（実装済み）**: 不純物バブルと同じ帰属分解（スピンチャネル + 軌道分解）を格子バブル χ₀^latt に適用。出力ファイル `DMFT_attrib_chi0_lattice.dat`。

8. **BSE chi 帰属（実装済み）**: BSE 解法後の全感受率 χ_full に対して帰属分解を適用。出力ファイル `DMFT_attrib_chi_full.dat`。不純物バブル / 格子バブル / BSE chi の3段階で帰属を比較することで、k 点和と頂点補正の効果を個別に同定できる。

9. **レベル間帰属比較（実装済み）**: 不純物バブル / 格子バブル / BSE chi_full の3段階のスピン/軌道帰属を一つのファイル `DMFT_attrib_comparison.dat` で比較する。各ボソン周波数 iΩₘ における:
   - IMP / LATT / BSE の値を並列出力
   - Δ_disp = LATT − IMP（k 点分散効果）
   - Δ_vtx = BSE − LATT（頂点補正効果）
   を全スピンチャネルおよび軌道ペアについて出力する。

10. **多原子格子バブル（実装済み）**: 全相関原子（lpawu == maxlpawu）の投影子からの寄与を累積した格子バブルに対して帰属分解を適用。

11. **k 点分解軌道ペア帰属（実装済み）**: Phase 17 の k 点スピンチャネルトレースを軌道ペアレベルまで拡張。各 k 点における S⁺S⁻ / S⁻S⁺ の軌道ペア (m, m') 分解を保持し、`DMFT_attrib_chi0_kpoint.dat` の Section 3-4 として出力。Section 4 では全 (k, m, m') の組み合わせを |χ₀^{+-}_{mm'}(k, 0)| の降順にソートしたグローバルランキング（上位20件）を提供。

12. **周波数依存帰属プロファイル（実装済み）**: 各ボソン Matsubara 周波数における支配的スピンチャネルと軌道ペアを自動同定する `DMFT_attrib_freqprofile_*.dat`。3つのセクション:
    - チャネル分率 vs 周波数（frac_sc, frac_pm, frac_mp の周波数依存性）
    - 支配的軌道ペア vs 周波数（各周波数で最大の |χ^{+-}_{mm'}| を持つ (m, m') ペア）
    - 周波数安定性診断（支配チャネルが全周波数で一貫するかの判定）

13. **光学伝導度スピンチャネル帰属（実装済み）**: バブル光学伝導度 Π_μν^bubble をスピン保存 / S⁺S⁻ / S⁻S⁺ の3チャネルに分解。KS バンド基底での電流行列要素 j_μ^{ab}(k) と格子 Green 関数 G^{ab}(k,iω) の4重積 Σ_{abcd} j_μ^{ab} G^{bc} j_ν^{cd} G^{da} の外殻添字 (a,d) のスピン帰属に基づく厳密な分類。出力ファイル `DMFT_optic_kernel.dat`。

14. **光学伝導度の軌道ペア分解（実装済み）**: chipsi 投影子を用いて速度行列要素と Green 関数を相関軌道空間に射影し、バブル光学伝導度 Π_μν を軌道ペア (m,m') ごとに分解。射影バブル Π^{proj} のカバー率、S⁺S⁻/S⁻S⁺ チャネルの軌道ペア分解、支配的 (m,m') ランキング、周波数依存の軌道キャラクター変化を出力。出力ファイル `DMFT_attrib_optic_orbital.dat`。

### 今後実装すべき帰属機能

11. ~~**k 点分解**: 逆格子空間での寄与の分布を出力する~~ → **Phase 17 で実装済み**
12. ~~**光学応答の帰属**: 光学伝導度 Π_μν を軌道対ごとに分解し、各吸収ピークの軌道起源を同定する（電流行列要素実装後）~~ → **Phase 20-21 でスピンチャネル帰属、Phase 24 で軌道ペア分解を実装済み。**
13. **原子間交差項の帰属**: 異なる原子間の交差寄与を格子バブルレベルで分解する（lpawu が異なる原子を含む系向け）
14. ~~**k 点分解の軌道分解**: 各 k 点における軌道ペア分解 S⁺S⁻ 感受率を出力する（Phase 17 の拡張）~~ → **Phase 18 で実装済み**
15. ~~**周波数依存帰属プロファイル**: 各ボソン周波数での支配チャネル/軌道ペアの同定~~ → **Phase 19 で実装済み**

### 帰属出力ファイル一覧

| ファイル名 | レベル | 内容 |
| --- | --- | --- |
| `DMFT_attrib_chi0_imp_atom{N}.dat` | 不純物（原子別） | 原子 N の chi0_imp のスピン/軌道帰属 |
| `DMFT_attrib_chi0_imp_total.dat` | 不純物（合算） | 全原子合算の chi0_imp のスピン/軌道帰属 |
| `DMFT_attrib_summary_imp.dat` | 不純物サマリー | 静的感受率・チャネル分率・収束診断・軌道ランキング |
| `DMFT_attrib_freqprofile_imp.dat` | 不純物プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_chi0_lattice.dat` | 格子バブル | 格子バブル chi0_latt のスピン/軌道帰属 |
| `DMFT_attrib_chi0_kpoint.dat` | 格子バブル（k分解） | k 点分解スピンチャネル + 軌道ペア分解 + ランキング |
| `DMFT_attrib_summary_latt.dat` | 格子サマリー | 格子レベルの静的感受率・軌道ランキング |
| `DMFT_attrib_freqprofile_latt.dat` | 格子プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_chi_full.dat` | BSE 全感受率 | BSE 補正後 chi_full のスピン/軌道帰属 |
| `DMFT_attrib_summary_bse.dat` | BSE サマリー | BSE レベルの静的感受率・軌道ランキング |
| `DMFT_attrib_freqprofile_bse.dat` | BSE プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_comparison.dat` | レベル間比較 | 不純物/格子/BSE の3段階比較（差分付き） |
| `DMFT_optic_kernel.dat` | バブル光学伝導度 | Π_μν テンソルのスピンチャネル分解（SC/S⁺S⁻/S⁻S⁺） |
| `DMFT_attrib_optic_orbital.dat` | 光学伝導度（軌道分解） | カバー率、軌道ペア分解、ランキング、周波数依存性 |

---

## 技術的リスクと未解決問題

1. **TRIQS 二粒子測定の統計誤差**: CT-HYB の二粒子量は一粒子量よりはるかにノイジーであり、Γ = χ₀⁻¹ − χ⁻¹ の行列反転が数値不安定になりうる。設計書に従い、平滑化やヒューリスティックな正則化は行わない。対処手段は測定回数増加・高精度演算・対称性投影に限定する。

2. **メモリコスト**: Mn 3d 全軌道（N_orb=5, N_σ=2）で N_ω=100 のフェルミ周波数を用いると、BSE 行列は各ボソン周波数で 10,000×10,000 の複素行列となる。N_Ω=50 で頂点全体 ~80 GB。分散メモリ並列化が不可避。

3. **実周波数 backend の不在**: これがない限り、厳密な意味での吸収スペクトル（α(ω)）は完成しない。Matsubara 軸応答は中間生成物として正当だが、「吸収スペクトル完成」とは呼べない。

4. **格子バブルの k 点並列化**: 現在の `compute_chi0_lattice` はシリアル実装であり、大規模 k 点メッシュでは計算時間が問題になる。O(nboson × niw_vertex × nkpt × norb_corr⁴ × mbandc²) の演算量であり、実用的な計算には MPI 分散が必須。

5. **chipsi 添字順序の検証**: `chipsi` と `matlu%mat` のスピノル軌道添字の順序が同一であることを前提としている。この前提は ABINIT の内部実装に基づくが、明示的な検証テストは未実施。実際の計算結果で `check_spinor_completeness` の出力を確認し、投影子の完全性を検証する必要がある。

6. **Green 関数の寿命問題**: 収束した Green 関数は `dmft_solve` のローカル変数であり、DMFT ループ終了時に破壊される。現在の実装では `dmft_solve` 内部で吸収計算を呼び出すことでこの問題を回避しているが、将来的には Green 関数の永続化（ファイル出力/読み込み）を検討すべきである。

7. **has_operks の可用性**: `compute_chi0_lattice` は `green_imp%oper(iw)%has_operks == 1` を前提とする。KS 基底データがないケース（例: 特定の DMFT ソルバーでの省メモリモード）では格子バブルがゼロになる。`ladd=.true.`（累積モード）の場合は既存の値を保持し、他の原子の寄与を壊さない設計とした。

8. **lpawu 不整合原子の除外**: 格子バブル計算では `lpawu /= maxlpawu` の原子を除外する。これは複合添字空間の次元が `maxlpawu` で定義されているためであり、異なる `lpawu` の原子の投影子を同じ空間で累積できないことに起因する。MnF₂（全 Mn 原子が lpawu=2）では問題にならないが、異種原子系（例: d + f 混合系）では一部の原子の寄与が失われる。この制約の解消には複合添字空間の拡張が必要。

---

## Phase 20-22: 電流頂点と光学バブル伝導度の帰属分解

### 19. 電流頂点モジュール（Phase 20: 完了）

**目的:** 光学伝導度バブル Π_μν^bubble の計算に必要な運動量行列要素 <ψ_a(k)|−i∇_μ|ψ_b(k)> を DMFT フレームワーク内で計算する

**背景:**

光学吸収スペクトルの帰属を行うには、感受率 χ のレベルだけでなく、光学伝導度 Π_μν のレベルでスピンチャネルと軌道ペアの分解が必要である。Π_μν の計算には電流（速度）行列要素 j_μ(k) が不可欠であり、これまでの実装ではこの部分がスケルトンのままであった。

設計書 Section 5.13-5.14 に基づき、電流演算子は
```
j_μ = −e ∂H/∂k_μ
```
であり、KS 固有状態基底では
```
<ψ_a(k)|j_μ|ψ_b(k)> = −e × <ψ_a(k)|−i∇_μ|ψ_b(k)>
```
で与えられる。PAW 形式では運動量行列要素は2つの寄与からなる：
1. **運動学的部分**: Σ_G (k+G)_μ × c*_a(G) × c_b(G) — 平面波部分
2. **PAW 増強部分**: −i × Σ_{ij,atom} cprj*_{a,i} × ∇_ij × cprj_{b,j} — PAW 球内補正

ここで ∇_ij = <φ_i|∇|φ_j> − <t̃φ_i|∇|t̃φ_j> は全電子と擬部分波の差分の行列要素であり、pawtab%nabla_ij に格納される。

**実装した内容:**

1. **`paw_dmft_type` への速度行列格納フィールドの追加**: `m_paw_dmft.F90` に以下を追加：
   - `psinablapsi_dmft(2, 3, mbandc, mbandc, nkpt, nsppol)` — 運動量行列要素の格納配列
   - `has_psinablapsi_dmft` — 計算状態フラグ（0: 未計算、1: 計算済み）
   - `destroy_dmft` に ABI_SFREE 追加

2. **`m_dmft_current_vertex.F90` モジュールの新規作成**: 以下の公開サブルーチンを実装：
   - `compute_psinablapsi_dmft`: cg（平面波係数）、cprj（PAW 投影子係数）、kg（G ベクトル指標）を用いて、相関バンド窓 [dmftbandi, dmftbandf] 内の全バンドペア (a,b) について運動量行列要素を計算

3. **運動学的部分の計算**: 各 k 点で (k+G)_μ ベクトルを構築し、平面波係数の共役積との内積を計算：
   ```
   <ψ_a|−i∇_μ|ψ_b> ← Σ_G (k+G)_μ × c*_a(G) × c_b(G)
   ```
   nspinor=2 の場合、同一 G ベクトルがスピン成分間で共有されるため、kpg_k 配列を 2*npw_k サイズで確保し、第2スピノル成分にも同じ G ベクトルをコピーする。

4. **PAW 増強部分の計算**: pawtab%nabla_ij が利用可能な場合（has_nabla フラグで判定）、`pawcprj_get` を用いて各 k 点の cprj データを抽出し、PAW 球内補正を追加：
   ```
   <ψ_a|−i∇_μ|ψ_b> += −i × Σ_{ij,atom} cprj*_{a,i} × nabla_ij(μ) × cprj_{b,j}
   ```
   nabla_ij が未計算の場合は運動学的部分のみを計算し、WARNING を出力する。これはヒューリスティックな省略ではなく、nabla_ij が `pawnabla_init` の呼び出しに依存するという技術的制約の反映である。

5. **`vtorho` からの呼び出し接続**: `m_vtorho.F90` に以下を追加：
   - `use m_dmft_current_vertex` の追加
   - `datafordmft` の呼び出し直後（cg, cprj, kg がメモリ上にある時点）で `compute_psinablapsi_dmft` を呼び出し
   - `dtset%dmft_resp_mode > 0` の条件で制御
   - `cryst_struc`（atindx1 用）、`mpi_enreg%comm_kpt`（MPI 通信用）、`mpi_enreg%proc_distrb`（プロセス分散用）を引数として渡す

**データフローの変更:**

```
vtorho
  ├── datafordmft → chipsi, eigen_dft を paw_dmft に格納
  ├── compute_psinablapsi_dmft → psinablapsi_dmft を paw_dmft に格納 [NEW]
  └── dmft_solve
        └── dmft_absorption_run → psinablapsi_dmft を利用して bubble conductivity を計算
```

**変更ファイル:**
- `src/65_paw/m_paw_dmft.F90` — psinablapsi_dmft フィールド追加、destroy 更新
- `src/68_dmft/m_dmft_current_vertex.F90` — 新規モジュール
- `src/79_seqpar_mpi/m_vtorho.F90` — use 文追加、compute_psinablapsi_dmft 呼び出し追加

### 20. バブル光学伝導度の実装（Phase 21: 完了）

**目的:** スケルトンだった `compute_bubble_conductivity` を完全実装し、スピンチャネル帰属分解を追加する

**実装した内容:**

1. **`optic_kernel_type` のスピンチャネル分解拡張**: 既存の `pi_bubble`, `pi_vertex`, `pi_total` に加え、以下を追加：
   - `pi_bubble_sc(nboson, ndir, ndir)` — スピン保存成分
   - `pi_bubble_pm(nboson, ndir, ndir)` — S⁺S⁻ スピン反転成分
   - `pi_bubble_mp(nboson, ndir, ndir)` — S⁻S⁺ スピン反転成分

2. **`compute_bubble_conductivity` の完全実装**: バブル光学伝導度を厳密に計算：
   ```
   Π_μν^bubble(iΩ_m) = −(1/βN_k) Σ_{k,n} Σ_{a,b,c,d} j_μ^{ab}(k) G^{bc}(k,iω_n) j_ν^{cd}(k) G^{da}(k,iω_n+iΩ_m)
   ```
   ここで：
   - j_μ^{ab}(k) は `paw_dmft%psinablapsi_dmft` から構築
   - G^{ab}(k,iω) は `green_imp%oper(iw)%ks(a,b,ikpt,isppol)` からアクセス
   - `has_psinablapsi_dmft /= 1` の場合はゼロを返し、WARNING を出力

3. **スピンチャネル分解の実装**: nspinor=2 の場合、4重ループの外殻添字 (a,d) のスピン帰属に基づいて各寄与をスピンチャネルに分類：
   - **スピン保存** (sc): a と d が同じスピンブロック（a ≤ N/2 かつ d ≤ N/2、または a > N/2 かつ d > N/2）
   - **S⁺S⁻** (pm): a がスピンアップブロック（a ≤ N/2）、d がスピンダウンブロック（d > N/2）
   - **S⁻S⁺** (mp): a がスピンダウン、d がスピンアップ

   この分類は `classify_spin_pair` 純粋関数で実装。

4. **`write_optic_kernel` の拡張**: 出力を3セクションに拡張：
   - **Section 1**: 全光学カーネル Π_μν（バブル + 全体）
   - **Section 2**: スピンチャネル分解（各ボソン周波数でのSC/PM/MP成分）
   - **Section 3**: iΩ=0 での対角成分のスピンチャネル分率

5. **ドライバの更新**: `dmft_absorption_run` の Stage 5 で `compute_bubble_conductivity` に `green_imp`, `niw_vertex`, `nspinor` を追加で渡す。

**物理的意味:**

バブル光学伝導度のスピンチャネル分解により、光学吸収スペクトルのどの部分がスピン反転遷移に起因するかを直接同定できる。具体的には：
- Π_μν^{S+S-} が大きい偏光方向では、その方向の光がスピン反転励起を効率的に誘起する
- スピン保存成分 Π_μν^{sc} と S⁺S⁻ 成分 Π_μν^{pm} の比率から、スピン反転遷移の相対的重要性を定量化できる
- テンソル成分（μ,ν）ごとの分解により、偏光方向に依存したスピン反転吸収の異方性を検出できる（MnF₂ の正方晶系で重要）

**出力ファイル:**

| ファイル名 | 内容 |
| --- | --- |
| `DMFT_optic_kernel.dat` | Matsubara 軸上の光学伝導度テンソル（スピンチャネル分解付き） |

**帰属の新機能:**

光学伝導度レベルの帰属により、感受率 χ のレベルの帰属（Phase 8-19 で実装済み）に加えて、光学吸収に直接対応する物理量でのスピン/軌道帰属が可能になった：

| 帰属レベル | 物理量 | 出力 |
| --- | --- | --- |
| 不純物バブル χ₀^imp | 局所感受率 | DMFT_attrib_chi0_imp_*.dat |
| 格子バブル χ₀^latt | 格子感受率 | DMFT_attrib_chi0_lattice.dat |
| BSE 全感受率 χ_full | 頂点補正込み感受率 | DMFT_attrib_chi_full.dat |
| **バブル光学伝導度 Π^bubble** | **電流−電流相関関数** | **DMFT_optic_kernel.dat** |

**変更ファイル:**
- `src/68_dmft/m_dmft_optic_kernel.F90` — compute_bubble_conductivity 完全実装、スピン分解追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — Stage 5 の引数更新
- `src/68_dmft/CMakeLists.txt` — m_dmft_current_vertex.F90 追加
- `src/68_dmft/abinit.src` — m_dmft_current_vertex.F90 追加

### 正直な到達点の評価（Phase 22 時点）

**現時点で完成しているもの:**
- 入力変数体系と整合性検査
- 全モジュールのデータ構造定義
- **不純物バブル χ₀^imp の完全な計算**
- **多原子サポート**
- **スピンチャネル帰属分解**（3チャネルへの厳密分解）
- **軌道分解帰属**
- **スピノル投影子の chipsi 接続**
- **格子バブル χ₀^latt の完全な計算**
- **多原子格子バブル**
- **格子レベル帰属分解**
- **レベル間帰属比較**
- **帰属サマリーの物理量抽出**
- **k 点分解格子バブル帰属**
- **k 点分解軌道ペア帰属**
- **周波数依存帰属プロファイル**
- 既約頂点抽出の行列演算
- 格子 BSE 解法の行列演算
- Matsubara 軸での出力フォーマット
- 高水準ドライバによる全ステージのオーケストレーション
- DMFT ループから吸収計算ドライバへの呼び出し接続
- **運動量行列要素 <ψ_a|−i∇|ψ_b> の計算**（運動学的部分 + PAW 増強部分）[NEW]
- **バブル光学伝導度 Π_μν^bubble の完全な計算**[NEW]
- **光学伝導度のスピンチャネル帰属分解**（SC/S⁺S⁻/S⁻S⁺）[NEW]

**現時点で完成していないもの:**
1. **局所二粒子相関関数 χ^imp の測定**: TRIQS/CT-HYB インターフェースへの連携口。**これが最大のボトルネック**。
2. **k 点並列化**: compute_bubble_conductivity、compute_psinablapsi_dmft、compute_optic_orbital_attrib のMPI並列化
3. **実周波数応答バックエンド**: dmft_resp_mode=2 の実装

---

## Phase 23-24: PAW nabla_ij 自動初期化と光学伝導度の軌道ペア分解

### 21. PAW nabla_ij の自動初期化（Phase 23: 完了）

**目的:** `dmft_resp_mode > 0` の場合に PAW 増強行列要素 nabla_ij を自動的に利用可能にする

**背景:**

Phase 20 で実装した `compute_psinablapsi_dmft` は、運動量行列要素 <ψ_a|−i∇|ψ_b> の計算において PAW 増強部分（∇_ij = <φ_i|∇|φ_j> − <t̃φ_i|∇|t̃φ_j>）を利用する。しかし、この PAW 増強行列要素は `pawnabla_init` によって計算され、これまでは optics 後処理（`m_paw_optics.F90`）、スクリーニングドライバ（`m_screening_driver.F90`）、BSE ドライバ（`m_bethe_salpeter.F90`）でのみ呼ばれていた。

DMFT 応答計算のフローでは、`vtorho` → `datafordmft` → `compute_psinablapsi_dmft` → `dmft_solve` → `dmft_absorption_run` の順序で実行されるため、`pawnabla_init` が事前に呼ばれていない場合、`pawtab%nabla_ij` が未計算のまま velocity 行列要素の PAW 補正が省略されてしまう。

これは機能の欠落であり、ユーザーに対して `prtnabla > 0` の設定を暗黙的に要求するものであった。

**実装した内容:**

1. **`m_vtorho.F90` への `use m_paw_onsite, only : pawnabla_init` の追加**: vtorho モジュールのスコープに `pawnabla_init` を導入。

2. **`compute_psinablapsi_dmft` 呼び出し前の自動初期化**: `dmft_resp_mode > 0` の条件ブロック内で、`pawtab(:)%has_nabla` が全原子型で 2 未満の場合に `pawnabla_init(psps%mpsang, ntypat, pawrad, pawtab)` を自動的に呼び出す。

3. **引数の取得**: `pawnabla_init` が必要とする引数はすべて `vtorho` のスコープ内で利用可能:
   - `psps%mpsang`: `pseudopotential_type` から取得（最大角運動量 + 1）
   - `ntypat`: `vtorho` の直接引数
   - `pawrad`: `vtorho` の引数 `pawrad(psps%ntypat*psps%usepaw)`
   - `pawtab`: `vtorho` の引数 `pawtab(ntypat*psps%usepaw)`

4. **冪等性の保証**: `has_nabla >= 2` のチェックにより、既に計算済みの場合は重複計算を避ける。`pawnabla_init` 自体も内部で既存の `nabla_ij` を解放してから再計算するため、二重呼び出しに対してもメモリリークは発生しない。

**変更ファイル:**
- `src/79_seqpar_mpi/m_vtorho.F90` — `use m_paw_onsite` 追加、`pawnabla_init` 条件呼び出し追加

### 22. 光学伝導度の軌道ペア分解（Phase 24: 完了）

**目的:** バブル光学伝導度 Π_μν を相関軌道ペア (m,m') ごとに分解し、吸収スペクトルの各寄与がどの d 軌道間遷移に帰属されるかを同定できるようにする

**背景:**

Phase 21 で実装したバブル光学伝導度のスピンチャネル分解（SC/S⁺S⁻/S⁻S⁺）は、スピン反転遷移の全体的な寄与を定量化するが、個々の軌道遷移（例: d_{xy} → d_{xz}）の寄与は識別できなかった。MnF₂ のような反強磁性体では、結晶場分裂により d 軌道間の遷移エネルギーが異なるため、吸収スペクトルの異なる部分が異なる軌道遷移に帰属される可能性がある。

この軌道ペア分解は、設計書 Section 5.14 の光学伝導度の定式化に基づき、chipsi 投影子を用いた相関軌道空間への射影によって実現する。

**実装した内容:**

1. **`optic_orbital_attrib_type` データ構造の追加**: `m_dmft_optic_kernel.F90` に新しいデータ型を追加:
   - `pi_orb_pm(nboson, ndir, ndir, ndim_orb, ndim_orb)` — S⁺S⁻ チャネルの軌道ペア (m,m') 分解
   - `pi_orb_mp(nboson, ndir, ndir, ndim_orb, ndim_orb)` — S⁻S⁺ チャネルの軌道ペア分解
   - `pi_projected_total(nboson, ndir, ndir)` — 投影されたバブル全体（カバー率計算用）
   - `pi_projected_sc(nboson, ndir, ndir)` — スピン保存成分

2. **`compute_optic_orbital_attrib` サブルーチンの追加**: chipsi 投影子を用いて速度行列要素と Green 関数を相関軌道空間に射影し、軌道ペア分解を計算:

   **射影の定式化:**
   ```
   J_μ^{αβ}(k) = Σ_{a,b} P_{α,a}(k) j_μ^{ab}(k) P*_{β,b}(k)  (射影速度行列)
   G^loc_{αβ}(k,iω) = Σ_{a,b} P_{α,a}(k) G^{ab}(k,iω) P*_{β,b}(k)  (局所 Green 関数)
   ```

   **射影バブルの軌道分解:**
   ```
   Π_μν^{proj,(α,δ)}(iΩ) = -(w_k/β) Σ_n Σ_{β,γ} J_μ^{αβ} G^loc_{βγ}(iω) J_ν^{γδ} G^loc_{δα}(iω+iΩ)
   ```

   ここで (α,δ) の外殻スピノル軌道添字がスピンチャネルと軌道ペアを決定する:
   - S⁺S⁻: α=(m,↑), δ=(m',↓) → 軌道ペア (m, m')
   - S⁻S⁺: α=(m,↓), δ=(m',↑) → 軌道ペア (m, m')

   **計算手順:**
   ```
   prod1 = J_μ × G^loc(iω)            [norb_corr × norb_corr]
   prod2 = prod1 × J_ν                 [norb_corr × norb_corr]
   contrib(α,δ) = prod2(α,δ) × G^loc(δ,α;iω+iΩ)  [要素ごとの積]
   ```

3. **`ladd` パラメータによる多原子サポート**: `compute_chi0_lattice` と同様に、最初の原子では `ladd=.false.`（ゼロ初期化して計算）、2番目以降の原子では `ladd=.true.`（既存値に累積）として動作。

4. **`write_optic_orbital_attrib` サブルーチンの追加**: 5つのセクションからなる出力ファイルを生成:

   - **Section 1: カバー率** — 投影された Π^proj_total と完全な KS 基底の Π_bubble の比。これにより、光学応答のうちどの程度が相関軌道遷移で説明されるかを直接定量化できる。カバー率が低い場合（例: < 0.5）、非相関バンドからの寄与が支配的であり、軌道分解の解釈に注意が必要であることを意味する。

   - **Section 2: S⁺S⁻ 軌道ペア分解（iΩ=0）** — 各方向 μ について、全軌道ペア (m,m') の Π^{S+S-}_{mm'} の実部・虚部・絶対値を出力。

   - **Section 3: S⁻S⁺ 軌道ペア分解（iΩ=0）** — 同様に S⁻S⁺ チャネル。

   - **Section 4: 支配的 S⁺S⁻ 軌道ペアランキング** — 方向平均した |Π^{S+S-}_{mm'}(0)| の降順にソートし、上位20件を出力。これにより、スピン反転吸収に最も大きく寄与する軌道遷移を直接同定できる。

   - **Section 5: S⁺S⁻ 周波数依存性** — 各ボソン周波数 iΩ_m における全軌道ペアの |Π^{S+S-}_{mm'}| を出力。xx 方向（μ=1）について出力する。軌道キャラクターがエネルギースケールに依存するかどうかを検出できる。

5. **ドライバ統合（Stage 5a）**: `dmft_absorption_run` 内の Stage 5（バブル光学伝導度計算）の直後に Stage 5a として挿入。多原子ループで全相関原子（lpawu == maxlpawu）に対して投影子を順次適用し、軌道帰属を累積する。

**物理的意味:**

この実装により、感受率 χ のレベルの帰属（Phase 8-19 で実装済み）に加えて、光学伝導度 Π_μν のレベルでも軌道ペア分解が可能になった。これは以下の物理的情報を提供する:

1. **吸収ピークの軌道起源**: 例えば、Π^{S+S-}_{d_{xy},d_{xz}} が Π^{S+S-}_{d_{z²},d_{x²-y²}} より大きい場合、スピン反転吸収の主要な寄与が t₂g 間の遷移に帰属される。

2. **偏光方向と軌道の対応**: テンソル成分 μ=x,y,z ごとの分解により、特定の偏光方向で特定の軌道遷移が活性化することを同定できる（選択則の反映）。

3. **カバー率による帰属の信頼性評価**: pi_projected_total / pi_bubble_total のカバー率が低い場合、光学応答が非相関バンドに支配されており、d 軌道軌道帰属の意味が限定的であることを客観的に判定できる。

4. **周波数依存の軌道キャラクター変化**: Section 5 のデータにより、異なるエネルギースケール（Matsubara 周波数）で支配的な軌道遷移が変化するかどうかを検出できる。

**注意点:**

- **射影による情報の損失**: 相関軌道空間への射影は非ユニタリ（P†P ≠ I）であるため、投影バブル Π^proj は完全な Π_bubble より小さい。カバー率はこの損失を定量化する。
- **軌道基底の選択依存性**: 軌道添字 m は PAW 投影子（chipsi）で定義される基底に依存する。球面調和関数基底では m=-l,...,+l に対応するが、これは結晶場固有状態（例: t₂g, e_g）とは一般に一致しない。結晶場固有基底への変換は将来の拡張。
- **これは Matsubara 軸上の解析であり、実周波数の吸収スペクトルの個々のピークに対する直接的な帰属ではない。**

**メモリ考慮:** `pi_orb_pm` の配列サイズは nboson × 3 × 3 × ndim_orb² の複素数。Mn 3d（ndim_orb=5）、nboson=50 で 50 × 9 × 25 × 16 bytes = 180 KB。無視できる。

**出力ファイル:**

| ファイル名 | 内容 |
| --- | --- |
| `DMFT_attrib_optic_orbital.dat` | 軌道ペア分解光学伝導度帰属（カバー率、S⁺S⁻/S⁻S⁺ 軌道分解、ランキング、周波数依存性） |

**変更ファイル:**
- `src/68_dmft/m_dmft_optic_kernel.F90` — `optic_orbital_attrib_type`、`init/destroy/compute/write_optic_orbital_attrib` 追加、`zgemm_wrapper`/`zgemm_ct_wrapper`/`zgemm_nn_small` ヘルパー追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — Stage 5a 追加、use 文更新、`optic_orb_attrib` ローカル変数追加

### 正直な到達点の評価（Phase 24 時点）

**現時点で完成しているもの:**
- 入力変数体系と整合性検査
- 全モジュールのデータ構造定義
- **不純物バブル χ₀^imp の完全な計算**
- **多原子サポート**
- **スピンチャネル帰属分解**（3チャネルへの厳密分解）
- **軌道分解帰属**
- **スピノル投影子の chipsi 接続**
- **格子バブル χ₀^latt の完全な計算**
- **多原子格子バブル**
- **格子レベル帰属分解**
- **レベル間帰属比較**
- **帰属サマリーの物理量抽出**
- **k 点分解格子バブル帰属**
- **k 点分解軌道ペア帰属**
- **周波数依存帰属プロファイル**
- 既約頂点抽出の行列演算
- 格子 BSE 解法の行列演算
- Matsubara 軸での出力フォーマット
- 高水準ドライバによる全ステージのオーケストレーション
- DMFT ループから吸収計算ドライバへの呼び出し接続
- **運動量行列要素 <ψ_a|−i∇|ψ_b> の計算**（運動学的部分 + PAW 増強部分）
- **バブル光学伝導度 Π_μν^bubble の完全な計算**
- **光学伝導度のスピンチャネル帰属分解**（SC/S⁺S⁻/S⁻S⁺）
- **PAW nabla_ij の自動初期化**（dmft_resp_mode > 0 で自動有効化）[NEW]
- **光学伝導度の軌道ペア分解**（chipsi 投影による相関軌道空間への射影と (m,m') 分解）[NEW]

**現時点で完成していないもの:**
1. **局所二粒子相関関数 χ^imp の測定**: TRIQS/CT-HYB インターフェースへの連携口。**これが最大のボトルネック**。これが解決しない限り、既約頂点 Γ はゼロのままであり、BSE の結果はバブル近似と等価である。
2. **k 点並列化**: `compute_bubble_conductivity`、`compute_psinablapsi_dmft`、`compute_optic_orbital_attrib` のMPI並列化。大規模 k 点メッシュでの実用計算には不可避。
3. **実周波数応答バックエンド**: `dmft_resp_mode=2` の実装。これがない限り、厳密な意味での吸収スペクトル α(ω) は完成しない。

### 次ステップ

#### 次ステップ 1: TRIQS 二粒子インターフェース（最優先・最大のボトルネック）

**目的:** TRIQS/CT-HYB の `measure_G2_iw_ph` から χ^imp を取得する

**必要な作業:**
1. `src/67_triqs_ext/triqs_cthyb_qmc.cpp` に `measure_G2_iw_ph` パラメータを渡す口を追加
2. TRIQS/CT-HYB の二粒子 Green 関数出力を Fortran 配列に変換するブリッジ関数（`ISO_C_BINDING` 使用）
3. 粒子正孔チャネルの Fourier 規約を設計書 Section 5.7 と厳密に一致させる
4. `m_forctqmc.F90` の `ctqmc_calltriqs` に二粒子測定のフラグと出力先を追加
5. 受信した χ^imp データを `chi_loc_type` に格納する変換ルーチン

**注意点:**
- TRIQS/CT-HYB の `measure_G2_iw_ph` は標準機能だが、ABINIT インターフェースがこの出力を受け取る口を持っていない
- C++ と Fortran の間の複素数配列の受け渡しは `ISO_C_BINDING` を用いる
- 二粒子量のメモリコストは一粒子量の O(N²_ω) 倍であり、並列化が不可避
- これが完成しない限り、Γ はゼロのままで BSE 結果はバブル近似と等価

#### 次ステップ 2: k 点並列化

**目的:** 格子バブル計算と光学伝導度計算のスケーラビリティ確保

**必要な作業:**
1. `compute_chi0_lattice` の k 点ループを MPI で分散
2. `compute_bubble_conductivity` の k 点ループを MPI で分散
3. `compute_optic_orbital_attrib` の k 点ループを MPI で分散
4. 既存の ABINIT MPI 分散機構（`mpi_enreg%my_kpttab` 等）を利用
5. k 点合算の MPI_ALLREDUCE 追加

#### 次ステップ 3: 実周波数応答バックエンド

**目的:** `dmft_resp_mode=2` の実装（最難関）

**候補手法:**
- 数値的解析接続は設計書で禁止されている（MaxEnt, Padé いずれも不可）
- 許されるのは:
  - (A) 実周波数の不純物応答を直接計算する補助ソルバー（例: NRG, ED, iPT）
  - (B) Lehmann 表示を明示的に用いる定式化
- 設計書の推奨は方針 A

**この段階に到達するまでは、Matsubara 軸応答（dmft_resp_mode=1）までの出力に留める。**

---

## 帰属出力ファイル一覧（Phase 24 時点）

| ファイル名 | レベル | 内容 |
| --- | --- | --- |
| `DMFT_attrib_chi0_imp_atom{N}.dat` | 不純物（原子別） | 原子 N の chi0_imp のスピン/軌道帰属 |
| `DMFT_attrib_chi0_imp_total.dat` | 不純物（合算） | 全原子合算の chi0_imp のスピン/軌道帰属 |
| `DMFT_attrib_summary_imp.dat` | 不純物サマリー | 静的感受率・チャネル分率・収束診断・軌道ランキング |
| `DMFT_attrib_freqprofile_imp.dat` | 不純物プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_chi0_lattice.dat` | 格子バブル | 格子バブル chi0_latt のスピン/軌道帰属 |
| `DMFT_attrib_chi0_kpoint.dat` | 格子バブル（k分解） | k 点分解スピンチャネル + 軌道ペア分解 + ランキング |
| `DMFT_attrib_summary_latt.dat` | 格子サマリー | 格子レベルの静的感受率・軌道ランキング |
| `DMFT_attrib_freqprofile_latt.dat` | 格子プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_chi_full.dat` | BSE 全感受率 | BSE 補正後 chi_full のスピン/軌道帰属 |
| `DMFT_attrib_summary_bse.dat` | BSE サマリー | BSE レベルの静的感受率・軌道ランキング |
| `DMFT_attrib_freqprofile_bse.dat` | BSE プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_comparison.dat` | レベル間比較 | 不純物/格子/BSE の3段階比較（差分付き） |
| `DMFT_optic_kernel.dat` | バブル光学伝導度 | Π_μν テンソルのスピンチャネル分解（SC/S⁺S⁻/S⁻S⁺） |
| `DMFT_attrib_optic_orbital.dat` | 光学伝導度（軌道分解）[NEW] | カバー率、軌道ペア分解、ランキング、周波数依存性 |

---

## Phase 25: TRIQS 二粒子測定インターフェース（完了）

### 23. TRIQS/CT-HYB G2_iw_ph 測定インターフェース（Phase 25: 完了）

**目的:** TRIQS/CT-HYB ソルバーの `measure_G2_iw_ph` 機能を通じて、不純物二粒子相関関数 χ^imp を取得するためのインターフェースを実装する。これにより、既約頂点 Γ がゼロではない物理的に意味のある値を持つようになり、BSE の結果がバブル近似を超えた多体効果を反映する。

**背景:**

Phase 24 までの実装で、バブル近似レベルの全計算パイプライン（不純物バブル → 格子バブル → BSE → 光学伝導度 → 帰属分解）は完成していた。しかし、BSE に入力される既約頂点 Γ^imp は常にゼロであり（χ^imp が未測定のため）、BSE の出力はバブル近似と等価であった。

設計書 Section 5.7 で定義される局所二粒子相関関数:

```
χ^imp_{αβγδ}(iω_n, iω_n'; iΩ_m) = -<T c†_α(τ₁) c_β(τ₂) c†_γ(τ₃) c_δ(0)>_connected
```

は、TRIQS/CT-HYB の `measure_G2_iw_ph` で測定される G2_iw_ph と符号を除いて同一である（χ = -G2）。この符号規約は Rohringer et al., Rev. Mod. Phys. 90, 025003 (2018) に準拠する。

**実装した内容:**

1. **`paw_dmft_type` への G2 格納フィールドの追加** (`m_paw_dmft.F90`):
   - `has_chi_imp_g2` — 測定状態フラグ（0: 未測定、1: 測定済み）
   - `chi_imp_g2_nboson` — 測定に使用したボソン Matsubara 周波数数
   - `chi_imp_g2_niw` — 測定に使用したフェルミ Matsubara 周波数数
   - `chi_imp_g2_norb` — フレーバー数（= nspinor × (2×lpawu+1)）
   - `chi_imp_g2_data(:)` — G2 データのフラット化配列

   データレイアウト:
   ```
   index = ((((iOm*niw + iw)*niw + iwp)*norb + α)*norb + β)*norb² + γ*norb + δ
   ```
   0-based インデックス。配列サイズ: nboson × niw² × norb⁴

2. **C++ ソルバーインターフェースの拡張** (`triqs_cthyb_qmc.cpp`, `triqs_cthyb_qmc.hpp`):
   - `ctqmc_triqs_run` 関数に5つの新パラメータを追加:
     - `bool measure_g2` — G2 測定の有効/無効
     - `int g2_n_bosonic` — ボソン周波数数
     - `int g2_n_fermionic` — フェルミ周波数数
     - `complex<double> *g2_data` — 出力データポインタ
     - `int g2_data_size` — 出力配列サイズ
   - ソルバーパラメータ設定:
     ```cpp
     paramCTQMC.measure_G2_iw_ph = true;
     paramCTQMC.measure_G2_n_bosonic = g2_n_bosonic;
     paramCTQMC.measure_G2_n_fermionic = g2_n_fermionic;
     ```
   - `solver.solve()` 後の G2 データ抽出: TRIQS の Block2Gf 構造からフラット配列へのブロック→フレーバー添字変換を実装。正の周波数のみを抽出:
     - ボソン: メッシュインデックス `g2_n_bosonic + iOm` (Ω=0 が `g2_n_bosonic`)
     - フェルミ: メッシュインデックス `g2_n_fermionic + iw` (ω₀ が `g2_n_fermionic`)

3. **Fortran ISO_C_BINDING の更新** (`triqs_interface_ctqmc.F90`):
   - `Ctqmc_triqs_run` の `bind(c)` インターフェースに新パラメータを追加:
     - `LOGICAL, VALUE :: measure_g2`
     - `INTEGER, VALUE :: g2_n_bosonic, g2_n_fermionic, g2_data_size`
     - `TYPE(C_PTR), VALUE :: g2_data`

4. **Fortran ラッパーの更新** (`m_forctqmc.F90`):
   - `ctqmc_calltriqs_c` に `optional :: measure_g2` フラグを追加
   - G2 測定パラメータを `paw_dmft` から読み取り
   - G2 バッファの動的確保と C ポインタへの変換
   - ソルバー呼び出し後、G2 データを `paw_dmft%chi_imp_g2_data` にコピー
   - `paw_dmft%has_chi_imp_g2 = 1` のフラグ設定

5. **G2 → chi_loc 変換ルーチン** (`m_dmft_two_particle.F90`):
   - `fill_chi_loc_from_g2` サブルーチンの新規追加
   - G2 フラット配列を chi_loc_type の複合添字形式に変換
   - 符号規約: χ = -G2（設計書 Section 5.7 の規約に準拠）
   - 次元整合性の検証（norb, niw, nboson の一致チェック）

6. **DMFT ドライバへの G2 測定呼び出し追加** (`m_dmft.F90`):
   - DMFT 自己無撞着ループ収束後、Green 関数が破壊される前に、G2 測定付きで不純物ソルバーを再実行:
     ```fortran
     if (dtset%dmft_resp_current_vertex == 1 .and. &
       & (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)) then
       paw_dmft%chi_imp_g2_nboson = dtset%dmft_resp_nboson
       paw_dmft%chi_imp_g2_niw = dtset%dmft_resp_niw_vertex
       call ctqmc_calltriqs_c(..., measure_g2=.true.)
     end if
     ```
   - 設計書 Section 9.2 に従い、収束した不純物浴を固定して二粒子測定のみを実行する

7. **吸収ドライバの Stage 2 更新** (`m_dmft_absorption_driver.F90`):
   - `paw_dmft%has_chi_imp_g2 == 1` の場合、G2 データから chi_loc を充填
   - 次元不整合（norb, nboson, niw のいずれかが一致しない場合）は WARNING を出力し、chi_loc をゼロのまま維持（ヒューリスティックな補間は行わない）
   - 測定された chi_imp を `DMFT_chi_imp_measured.dat` として出力
   - G2 データが利用できない場合は従来通り WARNING を出力し、バブル近似を使用

**データフローの変更:**

```
vtorho
  ├── datafordmft → chipsi, eigen_dft を paw_dmft に格納
  ├── compute_psinablapsi_dmft → 速度行列要素を paw_dmft に格納
  └── dmft_solve
        ├── DMFT self-consistent loop (既存)
        │     └── impurity_solve → ctqmc_calltriqs_c (通常の一粒子測定)
        ├── G2 measurement (NEW) → ctqmc_calltriqs_c(measure_g2=.true.)
        │     └── G2 data → paw_dmft%chi_imp_g2_data
        └── dmft_absorption_run
              └── Stage 2: fill_chi_loc_from_g2(chi_loc, paw_dmft%chi_imp_g2_data)
```

**Fourier 規約の整合性:**

TRIQS/CT-HYB の G2_iw_ph と設計書 Section 5.7 の χ^imp は同一の Fourier 規約を使用している:
- τ₄ = 0 固定の粒子正孔チャネル
- Fourier 因子: e^{-iω_nτ₁} e^{i(ω_n+Ω_m)τ₂} e^{-i(ω_n'+Ω_m)τ₃}
- 唯一の差異は全体の符号: χ = -G2

この規約の整合性は、Rohringer et al. (2018) の定義と TRIQS/CT-HYB のドキュメントの両方を参照して確認した。

**出力ファイル:**

| ファイル名 | 内容 |
| --- | --- |
| `DMFT_chi_imp_measured.dat` | TRIQS から測定された χ^imp（G2 データが利用可能な場合のみ出力） |

**変更ファイル:**
- `src/65_paw/m_paw_dmft.F90` — G2 格納フィールド追加、destroy 更新
- `src/67_triqs_ext/triqs_cthyb_qmc.cpp` — G2 測定パラメータ、データ抽出追加
- `src/67_triqs_ext/triqs_cthyb_qmc.hpp` — 関数宣言更新
- `src/67_triqs_ext/triqs_interface_ctqmc.F90` — ISO_C_BINDING インターフェース更新
- `src/68_dmft/m_forctqmc.F90` — measure_g2 オプション追加、G2 バッファ管理
- `src/68_dmft/m_dmft_two_particle.F90` — fill_chi_loc_from_g2 追加
- `src/68_dmft/m_dmft.F90` — G2 測定呼び出し追加
- `src/68_dmft/m_dmft_absorption_driver.F90` — Stage 2 更新、use 文更新

### 正直な到達点の評価（Phase 25 時点）

**現時点で完成しているもの:**
- 入力変数体系と整合性検査
- 全モジュールのデータ構造定義
- **不純物バブル χ₀^imp の完全な計算**
- **多原子サポート**
- **スピンチャネル帰属分解**（3チャネルへの厳密分解）
- **軌道分解帰属**
- **スピノル投影子の chipsi 接続**
- **格子バブル χ₀^latt の完全な計算**
- **多原子格子バブル**
- **格子レベル帰属分解**
- **レベル間帰属比較**
- **帰属サマリーの物理量抽出**
- **k 点分解格子バブル帰属**
- **k 点分解軌道ペア帰属**
- **周波数依存帰属プロファイル**
- 既約頂点抽出の行列演算
- 格子 BSE 解法の行列演算
- Matsubara 軸での出力フォーマット
- 高水準ドライバによる全ステージのオーケストレーション
- DMFT ループから吸収計算ドライバへの呼び出し接続
- **運動量行列要素 <ψ_a|−i∇|ψ_b> の計算**
- **バブル光学伝導度 Π_μν^bubble の完全な計算**
- **光学伝導度のスピンチャネル帰属分解**
- **PAW nabla_ij の自動初期化**
- **光学伝導度の軌道ペア分解**
- **TRIQS G2_iw_ph 測定インターフェース** [NEW]
- **G2 → chi_loc 変換** [NEW]
- **DMFT 収束後の G2 測定呼び出し** [NEW]

**インターフェースは完成しているが、実際の動作確認には以下が必要:**

1. **TRIQS/CT-HYB ライブラリのコンパイル環境**: G2 測定の C++ コードは `HAVE_TRIQS_INTERNAL || HAVE_TRIQS_v3_2` のプリプロセッサガードに囲まれており、TRIQS なしでもコンパイルは通るが、実際の G2 測定は TRIQS 環境でのみ動作する。

2. **TRIQS G2_iw_ph API の詳細検証**: TRIQS バージョン間で G2_iw_ph のブロック構造やメッシュ規約に差異がある可能性がある。特に:
   - `solver.G2_iw_ph` のアクセスパス（ポインタ vs 値返却）
   - Block2Gf のデータレイアウト（行優先 vs 列優先）
   - ボソン/フェルミメッシュのインデックス規約（0-based vs 1-based）
   これらは TRIQS ソースコードとの突き合わせが必要であり、本実装は TRIQS 3.x のドキュメントに基づいている。

3. **メモリ管理**: Mn 3d（norb=10）で nboson=10, niw=20 の場合、G2 データサイズは 10 × 20² × 10⁴ = 40M 複素数 = 640 MB。これは実用的な上限に近く、大規模なパラメータセットではメモリ不足になる可能性がある。

**現時点で完成していないもの:**

1. **TRIQS G2_iw_ph API の実テスト**: コンパイル環境の制約により、G2 データ抽出ロジックの実行時検証は未完了。特に Block2Gf のデータアクセスパス `(*solver.G2_iw_ph)(ib1, ib2).data()(...)` の正確性は TRIQS 環境でのテストが必要。

2. **k 点並列化**: `compute_bubble_conductivity`、`compute_chi0_lattice`、`compute_optic_orbital_attrib` の MPI 並列化。大規模 k 点メッシュでの実用計算には不可避。

3. **実周波数応答バックエンド**: `dmft_resp_mode=2` の実装。これがない限り、厳密な意味での吸収スペクトル α(ω) は完成しない。

---

## 次ステップで実装すべきこと

### 次ステップ 0: TRIQS G2 インターフェースの実テスト（最優先）

**目的:** Phase 25 で実装した G2 測定インターフェースを TRIQS 環境で実際にテストし、データ抽出が正しく動作することを確認する。

**具体的な検証項目:**

1. `solver.G2_iw_ph` へのアクセスが正しいこと（ポインタデリファレンスとブロックアクセス）
2. ボソン/フェルミメッシュのインデックスオフセットが正しいこと:
   - ボソン: メッシュインデックス `g2_n_bosonic` が Ω=0 に対応すること
   - フェルミ: メッシュインデックス `g2_n_fermionic` が ω₀ に対応すること
3. `data()` テンソルの添字順序が (bos, fer, fer, o1, o2, o3, o4) であること
4. `flavor_list` によるブロック→グローバルフレーバー添字変換が正しいこと
5. 符号規約 χ = -G2 の検証: 既知の原子極限解（U=0 で χ = χ₀）との比較
6. メモリ使用量の実測（小さなパラメータセットでの確認）

**必要条件:** TRIQS/CT-HYB がインストールされた計算環境と MnF₂ テスト入力ファイル

### 次ステップ 1: k 点並列化

**目的:** 格子バブル計算と光学伝導度計算のスケーラビリティ確保

**必要な作業:**
1. `compute_chi0_lattice` の k 点ループを MPI で分散
2. `compute_bubble_conductivity` の k 点ループを MPI で分散
3. `compute_optic_orbital_attrib` の k 点ループを MPI で分散
4. 既存の ABINIT MPI 分散機構（`mpi_enreg%my_kpttab` 等）を利用
5. k 点合算の MPI_ALLREDUCE 追加
6. k 点帰属データの MPI 集約

**技術的考慮:**
- ABINIT の既存 k 点並列化は `proc_distrb` テーブルで管理されている
- `compute_chi0_lattice` のフェルミ周波数ループは k 点ループの内側にあり、通信回数の最小化が必要
- 帰属データ（kpt_attrib）は各 k 点で独立に計算されるため、並列化は自然

### 次ステップ 2: 実周波数応答バックエンド（最難関）

**目的:** `dmft_resp_mode=2` の実装

**候補手法:**
- 数値的解析接続は設計書で禁止されている（MaxEnt, Padé いずれも不可）
- 許されるのは:
  - (A) 実周波数の不純物応答を直接計算する補助ソルバー（例: NRG, ED, iPT）
  - (B) Lehmann 表示を明示的に用いる定式化
- 設計書の推奨は方針 A

**この段階に到達するまでは、Matsubara 軸応答（dmft_resp_mode=1）までの出力に留める。**

---

## 帰属出力ファイル一覧（Phase 25 時点）

| ファイル名 | レベル | 内容 |
| --- | --- | --- |
| `DMFT_attrib_chi0_imp_atom{N}.dat` | 不純物（原子別） | 原子 N の chi0_imp のスピン/軌道帰属 |
| `DMFT_attrib_chi0_imp_total.dat` | 不純物（合算） | 全原子合算の chi0_imp のスピン/軌道帰属 |
| `DMFT_attrib_summary_imp.dat` | 不純物サマリー | 静的感受率・チャネル分率・収束診断・軌道ランキング |
| `DMFT_attrib_freqprofile_imp.dat` | 不純物プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_chi0_lattice.dat` | 格子バブル | 格子バブル chi0_latt のスピン/軌道帰属 |
| `DMFT_attrib_chi0_kpoint.dat` | 格子バブル（k分解） | k 点分解スピンチャネル + 軌道ペア分解 + ランキング |
| `DMFT_attrib_summary_latt.dat` | 格子サマリー | 格子レベルの静的感受率・軌道ランキング |
| `DMFT_attrib_freqprofile_latt.dat` | 格子プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_chi_full.dat` | BSE 全感受率 | BSE 補正後 chi_full のスピン/軌道帰属 |
| `DMFT_attrib_summary_bse.dat` | BSE サマリー | BSE レベルの静的感受率・軌道ランキング |
| `DMFT_attrib_freqprofile_bse.dat` | BSE プロファイル | 周波数依存帰属プロファイル（支配チャネル/軌道） |
| `DMFT_attrib_comparison.dat` | レベル間比較 | 不純物/格子/BSE の3段階比較（差分付き） |
| `DMFT_optic_kernel.dat` | バブル光学伝導度 | Π_μν テンソルのスピンチャネル分解（SC/S⁺S⁻/S⁻S⁺） |
| `DMFT_attrib_optic_orbital.dat` | 光学伝導度（軌道分解） | カバー率、軌道ペア分解、ランキング、周波数依存性 |
| `DMFT_chi_imp_measured.dat` | 測定 χ^imp [NEW] | TRIQS G2_iw_ph から変換された不純物二粒子相関関数 |
