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

### 正直な到達点の評価

**現時点で完成しているもの:**
- 入力変数体系と整合性検査
- 全モジュールのデータ構造定義
- 既約頂点抽出の行列演算（Γ = χ₀⁻¹ − χ⁻¹）
- 格子 BSE 解法の行列演算（χ = [χ₀⁻¹ − Γ]⁻¹）
- Matsubara 軸での出力フォーマット
- 高水準ドライバによる全ステージのオーケストレーション

**現時点で完成していないもの（スケルトンのままの部分）:**
1. **不純物 Green 関数の抽出**: `paw_dmft_type` 内の Green 関数データと `chi_loc_type` の接続
2. **局所二粒子相関関数の測定**: TRIQS/CT-HYB インターフェース（`triqs_cthyb_qmc.cpp`）への `measure_G2_iw_ph` 連携口
3. **格子 Green 関数の k 点ループ**: G(k, iωₙ) から格子バブル χ₀ˡᵃᵗᵗ への変換
4. **電流（速度）行列要素の計算**: PAW 補正を含むスピノル電流頂点 j_μ(k)
5. **実周波数応答 backend**: これは設計上の最難関であり、未実装

---

## 次ステップで実装すべきこと

### 次ステップ 1: 不純物 Green 関数の接続

**目的:** `compute_chi0_imp` の完全な実装

**必要な作業:**
1. `m_green.F90` の `green_type` から不純物 Green 関数 G^imp(iωₙ) を取得するインターフェースを作成
2. `green_type` の `oper_type` 配列（周波数ごとの行列）を chi_loc_type の複合添字空間にマッピング
3. chi₀^imp_{αβγδ}(iωₙ, iω_{n'}; iΩₘ) = −β δ_{nn'} G^imp_{δα}(iωₙ) G^imp_{βγ}(iωₙ + iΩₘ) の完全実装
4. 不純物バブルの自己無撞着性検証（χ → χ₀ at U→0）

**具体的な変更箇所:**
- `src/68_dmft/m_dmft_two_particle.F90` の `compute_chi0_imp` サブルーチン
- `m_green.F90` の `green_type` からのデータ抽出ユーティリティ追加

### 次ステップ 2: TRIQS 二粒子測定インターフェース

**目的:** CT-HYB ソルバーから χ^imp を取得する連携口

**必要な作業:**
1. `src/67_triqs_ext/triqs_cthyb_qmc.cpp` に `measure_G2_iw_ph` パラメータを渡す口を追加
2. TRIQS/CT-HYB の二粒子 Green 関数出力を Fortran 配列に変換するブリッジ関数
3. 粒子正孔チャネルの Fourier 規約を設計書 Section 5.7 と厳密に一致させる
4. `m_forctqmc.F90` の `ctqmc_calltriqs` に二粒子測定のフラグと出力先を追加

**注意点:**
- TRIQS/CT-HYB の `measure_G2_iw_ph` は標準機能だが、ABINIT インターフェースがこの出力を受け取る口を持っていない
- C++ と Fortran の間の複素数配列の受け渡しは `ISO_C_BINDING` を用いる
- 二粒子量のメモリコストは一粒子量の O(N²_ω) 倍であり、並列化が不可避

### 次ステップ 3: 格子 Green 関数と格子バブルの接続

**目的:** `compute_chi0_lattice` の完全な実装

**必要な作業:**
1. 各 k 点で格子 Green 関数 G(k, iωₙ) を構築（自己エネルギー埋め込み込み）
2. スピノル投影子を用いて相関部分空間へ射影
3. q=0 のバブル: χ₀(q=0, iωₙ, iΩₘ) = −(1/Nₖ) Σₖ G(k, iωₙ) G(k, iωₙ+iΩₘ)
4. k 点並列化（既存の MPI 分散機構を利用）

### 次ステップ 4: PAW 電流行列要素のスピノル拡張

**目的:** `compute_bubble_conductivity` の完全な実装

**必要な作業:**
1. `m_paw_optics.F90` の速度行列要素計算をスピノル基底に拡張
2. 非局所ポテンシャル補正 j^PAW_μ = j^local_μ + (i/ℏ)[V_NL, r_μ] のスピン混合成分を保持
3. テンソルの全3×3成分を保持（等方近似は行わない）

### 次ステップ 5: 実周波数応答 backend

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

### 現時点で組み込まれた帰属の仕組み

1. **スピンチャネル分離**: `dmft_resp_spinflip=1` の有無で計算を行い、差分を取ることで、スピン反転励起に起因する吸収強度を分離できる構造にしている

2. **バブル vs 頂点補正の分離**: `optic_kernel_type` に `pi_bubble`, `pi_vertex`, `pi_total` を別々に保存する設計とした。これにより、頂点補正（多体効果）の寄与を定量的に分離できる

3. **テンソル成分の完全保持**: 光学伝導度テンソル σ_μν の 3×3 成分すべてを保持する。MnF₂のルチル型正方晶では σ_xx = σ_yy ≠ σ_zz であり、偏光方向による吸収の異方性を直接出力する

### 今後実装すべき帰属機能

4. **軌道分解**: 格子バブルの k 合計の際に、特定の相関軌道ペア (α,β) からの寄与を分離出力する
5. **スピン-軌道分解**: 複合添字 α=(m,σ) を m と σ に分解し、各遷移の軌道的・スピン的性格を同定する
6. **k 点分解**: 逆格子空間での寄与の分布を出力する

---

## 技術的リスクと未解決問題

1. **TRIQS 二粒子測定の統計誤差**: CT-HYB の二粒子量は一粒子量よりはるかにノイジーであり、Γ = χ₀⁻¹ − χ⁻¹ の行列反転が数値不安定になりうる。設計書に従い、平滑化やヒューリスティックな正則化は行わない。対処手段は測定回数増加・高精度演算・対称性投影に限定する。

2. **メモリコスト**: Mn 3d 全軌道（N_orb=5, N_σ=2）で N_ω=100 のフェルミ周波数を用いると、BSE 行列は各ボソン周波数で 10,000×10,000 の複素行列となる。N_Ω=50 で頂点全体 ~80 GB。分散メモリ並列化が不可避。

3. **実周波数 backend の不在**: これがない限り、厳密な意味での吸収スペクトル（α(ω)）は完成しない。Matsubara 軸応答は中間生成物として正当だが、「吸収スペクトル完成」とは呼べない。
