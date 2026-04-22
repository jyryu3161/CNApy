# CNApy 변경 사항 (Changelog)

이 문서는 원본 CNApy [[cnapy-org/CNApy v1.2.7]](https://github.com/cnapy-org/CNApy)를 베이스로 이 포크에서 추가·변경된 모든 내용을 시간 역순으로 기록합니다.
변경 유형은 [Keep a Changelog](https://keepachangelog.com/) 규약에 따라 **Added / Changed / Fixed / Removed**로 분류했습니다.

---

## [Unreleased] — 2026-04-21

### Removed
- **GECKO dead dialog 파일 정리**: 통합 다이얼로그 도입 후 사용되지 않던 이전 다이얼로그 6개 삭제 (−1,774 lines)
  - `gecko_dialog.py`, `enzyme_usage_dialog.py`, `proteomics_dialog.py`, `kcat_whatif_dialog.py`, `unikp_dialog.py`, `unikp_predictor.py`

### Changed
- `ecmodel_data.py`의 stale docstring을 `gecko_unified_dialog` 기준으로 정정

---

## [0.3.0] GECKO 효소 제약 모델링 — 2026-04-20

### Added
- **GECKO 3.0 기반 ecModel 빌드**: Analysis 메뉴 > "Enzyme-Constrained Model (GECKO)…"
  - Full 및 Light 두 formulation 지원
  - GECKO 3.0 sign convention (v2 부호 규약) 정확히 구현
  - Isozyme 및 pseudoreaction 처리
  - Protein pool exchange 및 usage 반응 자동 생성
  - 파일: `cnapy/ecmodel/ecmodel_builder.py`, `expansion.py`, `ec_structure.py`, `ecmodel_data.py`, `exceptions.py`
- **통합 워크플로우 다이얼로그**: 사이드바 네비게이션으로 4개 페이지 통합
  - Build ecModel / Enzyme Usage Report / Proteomics Integration / kcat What-if Analysis
  - 파일: `cnapy/gui_elements/gecko_unified_dialog.py`
- **GECKO YAML 저장/로드**: `!!omap` 태그 기반 ordered mapping, plain GEM과 ecModel YAML 모두 지원, `.cna` 프로젝트 ZIP에 `ecmodel_data.json`으로 영속
  - v1→v2 sign convention 자동 마이그레이션 포함
  - 파일: `cnapy/ecmodel/yaml_io.py`
- **File > New project from YAML**: GECKO YAML을 직접 불러와 새 프로젝트 시작
- **File > New project from GECKO example**: Human-GEM / Yeast-GEM 번들 데이터셋 즉시 실행
  - Human-GEM: `human-GEM.yml`, DLKcat.tsv, uniprot.tsv
  - Yeast-GEM: `yeast-GEM.yml`, customKcats.tsv, uniprot.tsv
  - 파일: `cnapy/data/examples/gecko/`
- **GECKO 테스트 수트** (124개 테스트): Stage 1 signs/isozymes/pseudoreactions, YAML round-trip, CNA migration, paper parity, 번들 예제 빌드 검증
  - 파일: `cnapy/tests/test_ecmodel/`

### Fixed
- NaN kcat/MW 값이 stoichiometric computation에 들어가기 전 reject
- 외부 GECKO YAML 관용 처리 및 usage 반응 중복 추가 방지
- Isozyme / pseudoreaction 버그 수정 및 GECKO 3.0 sign convention 채택
- GECKO proteomics integration과 kcat constraint를 논문 공식과 일치하도록 조정

---

## [0.2.4] UI Refinement — 2026-03-20

### Changed
- **FilterableComboBox**를 Flux Response Analysis(FRA) 반응 선택기에 적용하여 substring 검색 지원

### Fixed
- 콤보박스/스핀박스 위에서 마우스 휠 스크롤로 값이 우발적으로 변경되는 문제 방지
  - 파일: `cnapy/utils.py` (`no_scroll` 헬퍼 추가)

---

## [0.2.3] Dynamic FBA 개선 — 2026-03-19

### Changed
- dFBA 기본 substrate 개수를 3개로 설정
- dFBA 다이얼로그에 **반응 선택기(picker)** 추가
- dFBA 테이블 UI 개선 및 **일괄 체크박스 토글** 기능 추가
- 파일: `cnapy/gui_elements/dynamic_fba_dialog.py`

---

## [0.2.2] Omics Gene Knockout — 2026-02-25

### Added
- **MOMA/ROOM 기준 플럭스(reference) 선택 다이얼로그**: 배치 분석 시 기준 플럭스를 선택하는 공용 UI
  - 파일: `cnapy/gui_elements/moma_room_reference_dialog.py`

---

## [0.2.1] Plot Customization — 2026-02-19

### Added
- **Plot Customization 다이얼로그**: FigureCanvasQTAgg 기반 플롯 공통 커스터마이제이션 UI
  - 제목, 축 라벨, 축 스케일(log/linear), 축 범위 조정
  - 9개 다이얼로그에 "Customize Plot" 버튼 통합:
    - Flux Response Analysis, FSEOF, FVSEOF, Gene Essentiality, Robustness, Flux Sampling, Flux Optimization, Yield Optimization, Phase plane/Yield space
  - 파일: `cnapy/gui_elements/plot_customization_dialog.py`

---

## [0.2.0] E-Flux2 & 의존성 정리 — 2026-02-13

### Added
- **E-Flux2 알고리즘** (Omics Integration에 LAD와 함께): True L2 norm (QP) 기반 플럭스 예측, pFBA/FBA fallback 지원
  - 메뉴 레이블: "Transcriptome-based Flux Prediction (LAD/E-Flux2)..."

### Fixed
- **GPR 규칙 재귀 평가**: 이전의 flat gene aggregation을 proper recursive tree 평가로 수정
  - OR → max, AND → min 트리 순회
- **누락된 의존성 복원**: 원본 CNApy에 있던 의존성을 `pyproject.toml`에 복원 (gurobipy 등)

---

## [0.1.2] Multi-condition Omics & UI — 2026-02-05

### Added
- **다중 조건 오믹스 통합 분석**: 여러 조건 파일을 동시에 로드하여 비교 분석
- **정렬 + Fold Change 컬럼** (오믹스 결과 테이블)
- **반응 선택 콤보박스에 실시간 텍스트 필터링** (substring 검색)
- **맵 → Reactions 탭 동기화**: 맵에서 반응을 클릭하면 Reactions 탭의 선택이 자동 동기화

### Changed
- 오믹스 통합 다이얼로그를 **non-modal**로 변경하여 다른 작업과 병행 가능

### Fixed
- LAD 분석의 optlang variable access 에러 수정

---

## [0.1.1] Strain Design 확장 — 2026-01-22

### Added
- **FSEOF Analysis** (Flux Scanning based on Enforced Objective Flux): 타겟 생산 플럭스와 상관된 반응을 식별하여 과발현/녹아웃 대사 엔지니어링 타겟 제안
  - 파일: `cnapy/gui_elements/fseof_dialog.py`
- **FVSEOF Analysis** (FVA 기반 FSEOF): 각 스캔 지점에서 FVA를 수행하여 보다 엄밀한 타겟 식별
  - 파일: `cnapy/gui_elements/fvseof_dialog.py`
- **Batch MOMA/ROOM Analysis**: 여러 knockout 시나리오에 대해 MOMA/ROOM 일괄 실행, 제약 조건 초기화 지원
  - 파일: `cnapy/gui_elements/batch_moma_room_dialog.py`
- **Gene Essentiality Analysis**: 유전자 필수성 체계적 스크리닝
  - 파일: `cnapy/gui_elements/gene_essentiality_dialog.py`
- **Robustness Analysis**: 파라미터 변동에 대한 모델 강건성 평가
  - 파일: `cnapy/gui_elements/robustness_analysis_dialog.py`
- **MOMA/ROOM 템플릿 플럭스 선택기**: 기준 플럭스 선택 UI
- **반응 목록 강화**: equation 컬럼 및 direction 컨트롤 추가
- **분석 다이얼로그 UX 개선**

### Fixed
- Batch MOMA/ROOM의 multiprocessing 버그
- 맵 배경 SVG 파일에서 quote 텍스트 제거
- 콘솔 → `print()` 참조 변경

### Removed
- Python 콘솔 및 관련 UI 정리 (향후 확장을 위한 리팩토링)

---

## [0.1.0] 초기 확장 — 2025-12

### Added

#### 분석 기능
- **ROOM (Regulatory On/Off Minimization)**: MILP solver 기반 유전자 녹아웃 후 플럭스 변화 최소화
  - 요구사항: CPLEX / Gurobi / GLPK 등 MILP solver
  - 파일: `cnapy/moma.py`
- **Linear MOMA**: 선형 MOMA 분석 및 외부 의존성의 선택적 동작
- **Flux Sampling**: GUI 다이얼로그 기반 플럭스 샘플링
  - 파일: `cnapy/flux_sampling.py`, `cnapy/gui_elements/flux_sampling_dialog.py`
- **Flux Response Analysis**: 타겟 반응 플럭스 스캔 + 제품 최대 생산률 플롯
  - 파일: `cnapy/gui_elements/flux_response_dialog.py`
- **Dynamic FBA (dFBA)**: FBA + ODE 커플링 시간 경과 시뮬레이션
  - 참고문헌: Mahadevan et al. 2002, Varma & Palsson 1994
  - 파일: `cnapy/gui_elements/dynamic_fba_dialog.py`
- **Omics Integration (LAD)**: Transcriptome 기반 Least Absolute Deviation 플럭스 예측
  - Gene expression 데이터 로드 (CSV/TSV/Excel)
  - Gene-to-reaction 매핑 (GPR rules)
  - 다양한 aggregation 방법 (min/max/mean/sum)
  - 파일: `cnapy/gui_elements/omics_integration_dialog.py`
- **설정 가능한 Auto Analysis 방법** (FBA / MOMA)

#### 모델 및 시나리오 관리
- **Model Management** 도구들:
  - GPR 정리 (중복 유전자 자동 탐지 및 정리)
  - Dead-end Metabolites 탐지
  - Blocked Reactions 탐지 (FVA 기반)
  - Orphan Reactions 탐지
  - Model Validation (질량/전하 균형, 바운드 오류 등)
  - 파일: `cnapy/gui_elements/model_management_dialog.py`
- **External Flux Data Loading**: CSV/TSV에서 reaction-flux 데이터 로드, 다중 조건 비교, Log2 Fold Change 히트맵(녹색=상향, 빨강=하향)
  - 파일: `cnapy/gui_elements/flux_data_dialog.py`
- **Scenario Templates & Bookmarks** (Ctrl+T): predefined 배양 조건 템플릿, 빠른 knockout 생성, 시나리오 북마크
  - 파일: `cnapy/gui_elements/scenario_templates_dialog.py`
- **Media Management** (Ctrl+M): 배양 배지 구성 관리
  - 파일: `cnapy/gui_elements/media_management_dialog.py`

#### AI 통합
- **LLM 기반 Strain Analysis**: ChatGPT / Google Gemini로 반응·유전자 존재 가능성 분석
  - OpenAI GPT-4o 등 및 Google Gemini Flash 지원
  - 웹 검색 기반 실시간 정보 활용
  - API 키 로컬 저장, 결과 JSON/CSV 내보내기
  - 파일: `cnapy/gui_elements/llm_analysis_dialog.py`

#### 지도 기능
- **PNG/SVG 이미지만으로 맵 생성**: JSON 파일 없이 이미지 파일만으로 CNApy 맵 생성
- **커스텀 반응 박스**: 모델에 없는 반응 ID도 맵에 박스로 추가하여 플럭스 값 표시
- 파일: `cnapy/gui_elements/central_widget.py`

#### UI/UX
- **Alt+Left-click / 컨텍스트 메뉴로 반응 Knockout 토글**
- **OptKnock 설명 개선** (Strain Design 다이얼로그):
  - Outer Objective 예시: `EX_succ_e` (숙신산 생산)
  - Inner Objective 예시: `BIOMASS` (성장)

### Fixed
- 예외 처리 개선 및 코드 품질 향상
- 여러 작은 버그 수정

---

## [Base] 원본 CNApy 1.2.7 — 2025-11-04

이 포크의 베이스. [`cnapy-org/CNApy` v1.2.7](https://github.com/cnapy-org/CNApy/releases/tag/v1.2.7).
원본에 포함된 기능은 이 Changelog에 기록하지 않으며, 원본 저장소의 릴리스 노트를 참고하세요.

---

## 라이센스

이 변경 사항들은 **Apache License 2.0** 하에 배포되며, 원본 CNApy 프로젝트의 일부입니다. 모든 변경 사항은 원본 CNApy 프로젝트의 라이센스와 호환됩니다.
