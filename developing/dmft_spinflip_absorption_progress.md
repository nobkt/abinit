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

| モジュール | ファイル | 主要な型/ルーチン | 状態 |
| --- | --- | --- | --- |
| `m_dmft_spinor_proj` | `m_dmft_spinor_proj.F90` | `spinor_proj_type`, `init_spinor_proj`, `check_spinor_completeness` | データ構造とインターフェース完了。chipsi からの投影子展開は接続待ち |
| `m_dmft_two_particle` | `m_dmft_two_particle.F90` | `chi_loc_type`, `compute_chi0_imp`, `write_chi_loc` | データ構造完了。不純物 Green 関数からの chi0 計算はスケルトン |
| `m_dmft_vertex` | `m_dmft_vertex.F90` | `vertex_irr_type`, `extract_vertex_irr`, `write_vertex_irr` | 完全な頂点抽出アルゴリズム実装済み（Γ = χ₀⁻¹ − χ⁻¹） |
| `m_dmft_lattice_bse` | `m_dmft_lattice_bse.F90` | `lattice_bse_type`, `compute_chi0_lattice`, `solve_lattice_bse` | BSE 解法アルゴリズム実装済み。格子バブルの k 点ループはスケルトン |
| `m_dmft_optic_kernel` | `m_dmft_optic_kernel.F90` | `optic_kernel_type`, `compute_bubble_conductivity`, `write_optic_kernel` | データ構造とI/O完了。電流行列要素の計算はスケルトン |
| `m_dmft_absorption_driver` | `m_dmft_absorption_driver.F90` | `dmft_absorption_run` | 全ステージのオーケストレーション完了 |

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

3. **出力**: 各ボソン Matsubara 周波数 iΩₘ における各チャネルの感受率トレースと軌道分解データを `DMFT_attrib_chi0_imp.dat` に出力。

**データ型 `spectral_attribution_type` のメンバ:**
- `chi_total(nboson)`: 全トレース
- `chi_spin_conserving(nboson)`: スピン保存成分
- `chi_spin_flip_pm(nboson)`: S⁺S⁻ 成分
- `chi_spin_flip_mp(nboson)`: S⁻S⁺ 成分
- `chi_pm_orbital(nboson, ndim_orb, ndim_orb)`: 軌道分解 S⁺S⁻
- `chi_mp_orbital(nboson, ndim_orb, ndim_orb)`: 軌道分解 S⁻S⁺

**スピノル軌道添字の規約:**
- nspinor=2 の場合: α = 1..ndim_orb は spin-up（軌道 m = α）、α = ndim_orb+1..2×ndim_orb は spin-down（軌道 m = α − ndim_orb）
- ここで ndim_orb = 2×lpawu + 1

**注意事項:**
- nspinor=2 を要求する（スピン反転分解は nspinor=1 では物理的に無意味）
- ヒューリスティックな分類やしきい値判定は一切行っていない
- 分解は数学的に厳密であり、chi_total = chi_spin_conserving + chi_spin_flip_pm + chi_spin_flip_mp が成立する

**変更ファイル:**
- `src/68_dmft/m_dmft_spectral_attribution.F90` — 新規作成

#### 7. 吸収計算ドライバの Green 関数接続と帰属統合（Phase 5: 完了）

**実装した内容:**

1. **ドライバインターフェースの変更**: `dmft_absorption_run` のシグネチャを `(dtset, paw_dmft, cryst_struc)` から `(dtset, paw_dmft, cryst_struc, green_imp)` に変更し、収束した Green 関数を直接受け取るようにした。

2. **帰属処理の統合**: Stage 1 の chi0 計算直後に Stage 1b として帰属分解を追加。`dmft_resp_spinflip=1` かつ `nspinor=2` の場合のみ実行。

3. **ファイル出力の追加**:
   - `DMFT_chi0_imp.dat` — 不純物バブルの複合添字行列要素
   - `DMFT_attrib_chi0_imp.dat` — スピンチャネル別・軌道別の帰属分解結果
   - `DMFT_optic_kernel.dat` — Matsubara 軸光学カーネル

4. **DMFT ループからの呼び出し接続**: `dmft_solve` に `dtset` 引数を追加し、DMFT ループ収束後・Green 関数破壊前に `dmft_absorption_run` を呼び出すコードを挿入した。

**変更ファイル:**
- `src/68_dmft/m_dmft_absorption_driver.F90` — インターフェース変更と帰属統合
- `src/68_dmft/m_dmft.F90` — `dtset` 引数追加と吸収ドライバ呼び出し
- `src/79_seqpar_mpi/m_vtorho.F90` — `dmft_solve` 呼び出しの引数更新

#### 8. ビルドシステム更新（Phase 5b: 完了）

- `src/68_dmft/CMakeLists.txt` — `m_dmft_spectral_attribution.F90` をアルファベット順で追加
- `src/68_dmft/abinit.src` — 同上

### 正直な到達点の評価

**現時点で完成しているもの:**
- 入力変数体系と整合性検査
- 全モジュールのデータ構造定義
- **不純物バブル χ₀^imp の完全な計算**（Green 関数から直接構築）
- **スピンチャネル帰属分解**（スピン保存 / S⁺S⁻ / S⁻S⁺ の3チャネルへの厳密分解）
- **軌道分解帰属**（各 d 軌道ペアからの S⁺S⁻ 寄与の個別出力）
- 既約頂点抽出の行列演算（Γ = χ₀⁻¹ − χ⁻¹）
- 格子 BSE 解法の行列演算（χ = [χ₀⁻¹ − Γ]⁻¹）
- Matsubara 軸での出力フォーマット
- 高水準ドライバによる全ステージのオーケストレーション
- **DMFT ループから吸収計算ドライバへの呼び出し接続**

**現時点で完成していないもの（スケルトンのままの部分）:**
1. **局所二粒子相関関数 χ^imp の測定**: TRIQS/CT-HYB インターフェース（`triqs_cthyb_qmc.cpp`）への `measure_G2_iw_ph` 連携口。これが完成しない限り、Γ は自明な値（ゼロ）のままであり、BSE の結果はバブル近似と等価。
2. **格子 Green 関数の k 点ループ**: G(k, iωₙ) から格子バブル χ₀ˡᵃᵗᵗ への変換
3. **電流（速度）行列要素の計算**: PAW 補正を含むスピノル電流頂点 j_μ(k)
4. **実周波数応答 backend**: これは設計上の最難関であり、未実装
5. **スピノル投影子 chipsi の接続**: `m_dmft_spinor_proj.F90` の投影子配列を `paw_dmft%chipsi` から取得する部分

---

## 次ステップで実装すべきこと

### 次ステップ 1: TRIQS 二粒子測定インターフェース（最優先）

**目的:** CT-HYB ソルバーから χ^imp を取得する連携口の作成

**現状の問題:** χ^imp がゼロのままであるため、既約頂点 Γ も全出力もバブル近似と区別がつかない。これが解決しない限り、DMFT 応答計算の本質的な部分（多体効果の取り込み）が機能しない。

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

### 次ステップ 2: 格子 Green 関数と格子バブルの接続

**目的:** `compute_chi0_lattice` の完全な実装

**必要な作業:**
1. 各 k 点で格子 Green 関数 G(k, iωₙ) を構築（自己エネルギー埋め込み込み）
2. スピノル投影子を用いて相関部分空間へ射影
3. q=0 のバブル: χ₀(q=0, iωₙ, iΩₘ) = −(1/Nₖ) Σₖ G(k, iωₙ) G(k, iωₙ+iΩₘ)
4. k 点並列化（既存の MPI 分散機構を利用）
5. **格子バブルの帰属分解**: 既に実装した `m_dmft_spectral_attribution` の枠組みを格子バブルにも適用する

### 次ステップ 3: PAW 電流行列要素のスピノル拡張

**目的:** `compute_bubble_conductivity` の完全な実装

**必要な作業:**
1. `m_paw_optics.F90` の速度行列要素計算をスピノル基底に拡張
2. 非局所ポテンシャル補正 j^PAW_μ = j^local_μ + (i/ℏ)[V_NL, r_μ] のスピン混合成分を保持
3. テンソルの全3×3成分を保持（等方近似は行わない）

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

### 今後実装すべき帰属機能

6. **格子バブルの帰属分解**: 不純物バブルと同じ帰属分解を格子バブル χ₀ˡᵃᵗᵗ に適用する（格子バブル接続後）
7. **k 点分解**: 逆格子空間での寄与の分布を出力する
8. **光学応答の帰属**: 光学伝導度 Π_μν を軌道対ごとに分解し、各吸収ピークの軌道起源を同定する

---

## 技術的リスクと未解決問題

1. **TRIQS 二粒子測定の統計誤差**: CT-HYB の二粒子量は一粒子量よりはるかにノイジーであり、Γ = χ₀⁻¹ − χ⁻¹ の行列反転が数値不安定になりうる。設計書に従い、平滑化やヒューリスティックな正則化は行わない。対処手段は測定回数増加・高精度演算・対称性投影に限定する。

2. **メモリコスト**: Mn 3d 全軌道（N_orb=5, N_σ=2）で N_ω=100 のフェルミ周波数を用いると、BSE 行列は各ボソン周波数で 10,000×10,000 の複素行列となる。N_Ω=50 で頂点全体 ~80 GB。分散メモリ並列化が不可避。

3. **実周波数 backend の不在**: これがない限り、厳密な意味での吸収スペクトル（α(ω)）は完成しない。Matsubara 軸応答は中間生成物として正当だが、「吸収スペクトル完成」とは呼べない。

4. **多原子系への拡張**: 現在の compute_chi0_imp は最初の相関原子のみを処理する。MnF₂のような複数の Mn サイトを持つ系では、全相関原子の寄与を扱う拡張が必要。ただし、DMFT では各不純物問題は独立に解かれるため、原子ごとの処理をループ化すれば対応可能。

5. **Green 関数の寿命問題**: 収束した Green 関数は `dmft_solve` のローカル変数であり、DMFT ループ終了時に破壊される。現在の実装では `dmft_solve` 内部で吸収計算を呼び出すことでこの問題を回避しているが、将来的には Green 関数の永続化（ファイル出力/読み込み）を検討すべきである。
