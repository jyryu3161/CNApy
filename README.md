# CNApy: 대사 모델링을 위한 통합 환경

[![Latest stable release](https://flat.badgen.net/github/release/cnapy-org/cnapy/stable)](https://github.com/cnapy-org/CNApy/releases/latest)
[![Last commit](https://flat.badgen.net/github/last-commit/cnapy-org/cnapy)](https://github.com/cnapy-org/CNApy/commits/master)
[![Open issues](https://flat.badgen.net/github/open-issues/cnapy-org/cnapy)](https://github.com/cnapy-org/CNApy/issues)
[![Gitter chat](https://flat.badgen.net/gitter/members/cnapy-org/community)](https://gitter.im/cnapy-org/community)

![CNApy screenshot](CNApy-intro.png)

## 소개

안녕하세요! CNApy [[논문]](https://doi.org/10.1093/bioinformatics/btab828)에 오신 것을 환영합니다. CNApy는 파이썬 기반의 그래픽 사용자 인터페이스(GUI)로, 대사 모델링을 쉽고 편리하게 수행할 수 있도록 도와주는 도구입니다.

CNApy를 사용하면 다음과 같은 작업들을 할 수 있어요:
*   **다양한 COBRA 분석**: 화학양론적 대사 모델을 이용한 제약 조건 기반 재구성 및 분석(COBRA)의 여러 일반적인 방법들을 수행할 수 있습니다.
*   **인터랙티브한 시각화**: COBRA 계산 결과를 *인터랙티브하고 편집 가능한* 대사 지도로 시각화할 수 있습니다. Escher 지도 [[GitHub]](https://escher.github.io/#/)[[논문]](<https://doi.org/10.1371/journal.pcbi.1004321>)도 지원해요!
*   **모델 생성 및 편집**: 반응, 대사산물, 유전자 등을 포함한 대사 모델을 직접 만들고 수정할 수 있습니다.
*   **효소 제약 모델링 (GECKO)**: GECKO 3.0 기반 효소 제약 모델(ecModel)을 빌드하고 단백질체 데이터를 통합할 수 있습니다.

모델을 불러오거나 내보낼 때는 널리 사용되는 SBML 표준 형식 [[사이트]](https://sbml.org/)[[논문]](https://www.embopress.org/doi/abs/10.15252/msb.20199110)을 지원합니다.

지원하는 COBRA 방법들(일부는 cobrapy [[GitHub]](https://github.com/opencobra/cobrapy)[[논문]](https://doi.org/10.1186/1752-0509-7-74)에서 제공)은 다음과 같습니다:

- **Flux Balance Analysis (FBA)** [[리뷰]](https://doi.org/10.1038/nbt.1614)
- **Flux Variability Analysis (FVA)** [[논문]](https://doi.org/10.1016/j.ymben.2003.09.002)
- **Yield optimization** (선형 분수 프로그래밍 기반) [[논문]](https://doi.org/10.1016/j.ymben.2018.02.001)
- **Phase plane analyses** (플럭스 및/또는 수율 최적화 포함 가능)
- **Flux Sampling**: 모델의 가능한 플럭스 분포를 샘플링하여 분석할 수 있습니다.
- **Linear MOMA**: 기준 플럭스 분포와의 차이를 최소화하는 선형 MOMA 분석을 수행할 수 있습니다.
- **ROOM (Regulatory On/Off Minimization)**: MILP solver를 사용하여 유전자 녹아웃 후 유의미한 플럭스 변화를 최소화하는 분석을 수행할 수 있습니다 (MILP solver 필요).
- **Dynamic FBA (dFBA)**: ODE 커플링을 통해 시간 경과에 따른 바이오매스와 외부 대사체 농도 변화를 시뮬레이션합니다.
- **Flux Response Analysis**: 타겟 반응의 플럭스를 스캔하면서 제품 반응의 최대 생산률을 플롯하는 분석을 수행할 수 있습니다.
- **FSEOF / FVSEOF**: 타겟 생산 플럭스와 상관된 반응을 자동으로 식별하여 대사 엔지니어링 타겟(과발현/녹아웃)을 제안합니다.
- **Gene Essentiality Analysis**: 유전자 필수성을 체계적으로 스크리닝합니다.
- **Robustness Analysis**: 파라미터 변동에 대한 모델 강건성을 평가합니다.
- **Batch MOMA/ROOM Analysis**: 다중 knockout 시나리오를 한 번에 배치 실행합니다.
- **Omics Integration (LAD / E-Flux2)**: Transcriptome 데이터를 기반으로 LAD (Least Absolute Deviation) 또는 E-Flux2 (True L2 norm QP) 방법을 사용하여 플럭스 분포를 예측하고, 다중 조건을 비교할 수 있습니다.
- **GECKO 효소 제약 모델링** [[GitHub]](https://github.com/SysBioChalmers/GECKO)[[논문]](https://www.nature.com/articles/s41596-023-00931-7): GECKO 3.0 논문 방법론에 따라 기존 GEM을 ecModel로 확장하고, 단백질체 데이터 통합 및 kcat What-if 분석을 수행할 수 있습니다.
- **실제 측정된 *in vivo* 플럭스 시나리오 적용**: 화학양론적으로 타당하게 만들며, 선택적으로 바이오매스 반응을 수정할 수도 있습니다 [[논문]](https://academic.oup.com/bioinformatics/article/39/10/btad600/7284109).
- **Elementary Flux Modes (EFM)** [[리뷰]](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/biot.201200269)
- **OptMDFpathway 기반 열역학적 방법** [[논문]](https://doi.org/10.1371/journal.pcbi.1006492)
- **고급 균주 설계 알고리즘**: StrainDesign [[GitHub]](https://github.com/klamt-lab/straindesign)[[논문]](https://doi.org/10.1093/bioinformatics/btac632) 통합을 통해 OptKnock [[논문]](https://doi.org/10.1002/bit.10803), RobustKnock [[논문]](https://doi.org/10.1093/bioinformatics/btp704), OptCouple [[논문]](https://doi.org/10.1016/j.mec.2019.e00087), 고급 Minimal Cut Sets [[논문]](https://doi.org/10.1371/journal.pcbi.1008110) 등을 지원합니다.

**→ CNApy 설치 방법이 궁금하시다면, [설치 옵션](#설치-옵션) 섹션을 확인해주세요.**

**→ CNApy의 다양한 기능에 대해 더 알고 싶으시다면, [문서 및 튜토리얼](#문서-및-튜토리얼) 섹션을 참고해주세요.**

**→ 질문이나 제안, 버그 신고는 언제든 환영합니다! [CNApy GitHub 이슈](https://github.com/cnapy-org/CNApy/issues), [CNApy GitHub 토론](https://github.com/cnapy-org/CNApy/discussions) 또는 [CNApy Gitter 채팅방](https://gitter.im/cnapy-org/community)을 이용해주세요.**

**→ CNApy를 인용하고 싶으시다면, [CNApy 인용 방법](#cnapy-인용-방법) 섹션을 봐주세요.**

**→ 개발자로서 CNApy에 기여하고 싶으시다면, [CNApy 개발에 기여하기](#cnapy-개발에-기여하기) 섹션을 확인해주세요.**

*참고*: MATLAB 기반의 유명한 *CellNetAnalyzer* (CNA)를 사용하고 싶으시다면(CNApy와는 호환되지 않아요), [CNA 웹사이트](https://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html)에서 다운로드하실 수 있습니다.

## 설치 옵션

CNApy를 설치하는 4가지 방법이 있습니다:

1.  **가장 쉬운 방법**: 윈도우, 리눅스, 맥OS용 설치 프로그램을 다운로드하여 설치하는 것입니다. 자세한 내용은 [CNApy 설치 프로그램 사용하기](#cnapy-설치-프로그램-사용하기)를 참고하세요.
2.  **파이썬 사용자라면**: 시스템에 Python 3.10(다른 버전은 안 돼요!)이 설치되어 있다면, 콘솔에서 `pip install cnapy`를 입력하여 간단히 설치할 수 있습니다. 설치 후에는 `cnapy` 또는 `python -m cnapy` 명령어로 실행할 수 있습니다.
3.  **conda/mamba 사용자라면**: `cnapy-1.2.7`이라는 환경을 만들어 설치할 수 있습니다.
    1.  `conda create --name cnapy-1.2.7 python=3.10 pip openjdk -c conda-forge` 실행
    2.  `conda activate cnapy-1.2.7` 실행
    3.  `pip install cnapy` 실행
    이후 `cnapy` 또는 `python -m cnapy`로 실행하세요. (참고: [cnapy conda 패키지](https://anaconda.org/cnapy/cnapy)는 라이선스 문제로 현재 업데이트되지 않고 있습니다.)
4.  **개발자라면**: git과 conda/mamba를 사용하여 저장소를 복제하고 설정하는 방법은 [CNApy 개발 환경 설정하기](#cnapy-개발-환경-설정하기) 섹션을 참고해주세요.

## 문서 및 튜토리얼

-   [CNApy 가이드](https://cnapy-org.github.io/CNApy-guide/)에서 주요 기능에 대한 정보를 확인하실 수 있습니다.
-   [CNApy 유튜브 채널](https://www.youtube.com/channel/UCRIXSdzs5WnBE3_uukuNMlg)에서 사용 영상을 보실 수 있습니다.
-   [CNApy 예제 프로젝트](https://github.com/cnapy-org/CNApy-projects/releases/latest)를 제공합니다. 가장 일반적인 *E. coli* 모델들이 포함되어 있으며, CNApy를 처음 시작할 때나 파일 메뉴를 통해 다운로드할 수 있습니다.
-   **GECKO 예제**: File > "New project from GECKO example" 메뉴에서 Human-GEM과 Yeast-GEM 번들 데이터셋으로 GECKO 워크플로우를 즉시 체험해볼 수 있습니다.

## CNApy 설치 프로그램 사용하기

이 설치 프로그램을 사용하면 윈도우, 리눅스, 맥OS에서 로컬 설치를 쉽게 할 수 있습니다.

*윈도우 사용자:*

-   [여기서](https://github.com/cnapy-org/CNApy/releases/download/v1.2.7/install_cnapy_here.bat) 윈도우용 설치 프로그램을 다운로드하세요.
-   CNApy를 설치하고 싶은 폴더에 파일을 넣으세요.
-   파일을 더블 클릭하여 설치를 진행하세요.
-   설치가 끝나면 바탕화면의 CNApy 아이콘이나 `cnapy-1.2.7` 폴더 안의 "RUN_CNApy.bat"을 더블 클릭하여 실행할 수 있습니다.

*리눅스 또는 맥OS 사용자:*

-   [여기서](https://github.com/cnapy-org/CNApy/releases/download/v1.2.7/install_cnapy_here.sh) 리눅스 & 맥OS용 설치 프로그램을 다운로드하세요.
-   CNApy를 설치하고 싶은 폴더에 파일을 넣으세요.
-   콘솔을 열고 해당 폴더로 이동한 뒤 `chmod u+x ./install_cnapy_here.sh`를 실행하여 실행 권한을 부여하세요. (또는 파일 속성에서 실행 가능으로 설정하세요.)
-   콘솔에서 `./install_cnapy_here.sh`를 실행하거나 파일을 더블 클릭하여 설치하세요.
-   설치가 완료되면 콘솔에서 `./run_cnapy.sh`를 실행하거나(해당 폴더로 이동 필요), `cnapy-1.2.7` 폴더 안의 "run_cnapy.sh"를 더블 클릭하여 실행할 수 있습니다.

기술적 참고: CNApy 설치 프로그램은 [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)를 사용합니다.

## CNApy 개발 환경 설정하기

*참고:* 이 섹션은 CNApy 개발에 참여하고 싶은 프로그래머 분들을 위한 것입니다. 일반 사용자라면 [설치 옵션](#설치-옵션)을 따라주세요.

CNApy 개발에 참여해주시는 모든 분들을 환영합니다! [기여 가이드](https://github.com/cnapy-org/CNApy/blob/master/CONTRIBUTING.md)에서 일반적인 지침을 확인해주세요. 기여해주신 모든 내용은 Apache 2.0 라이선스 하에 배포됩니다.

개발을 위해서는 uv [[GitHub]](https://github.com/astral-sh/uv)를 사용하여 의존성과 파이썬 버전을 관리하는 것을 추천합니다. conda/mamba를 사용할 수도 있지만, 의존성을 수동으로 설치해야 할 수 있습니다.

### uv 사용법

1.  uv가 설치되어 있는지 확인하세요. (pip, pipx 등으로 설치 가능)

```sh
# 예시
pip install uv # 또는
pipx install uv
```

2.  git을 사용하여 최신 cnapy 개발 버전을 가져오세요.

```sh
git clone https://github.com/cnapy-org/CNApy.git
```

3.  소스 디렉토리로 이동하여 CNApy를 실행하세요.

```sh
cd CNApy
uv run cnapy.py
```

uv가 자동으로 올바른 파이썬 버전과 의존성을 설치해줍니다. 만약 Java/JDK/JVM/jpype 관련 오류가 발생한다면, 시스템에 OpenJDK [[사이트]](https://openjdk.org/install/)를 설치해보세요.

## CNApy 인용 방법

연구에 CNApy를 사용하셨다면, 다음 논문을 인용해주시면 감사하겠습니다:

Thiele et al. (2022). CNApy: a CellNetAnalyzer GUI in Python for analyzing and designing metabolic networks.
*Bioinformatics* 38, 1467-1469, [doi.org/10.1093/bioinformatics/btab828](https://doi.org/10.1093/bioinformatics/btab828).

## 최근 변경 사항

이 버전에 새로 추가된 주요 기능은 **GECKO 효소 제약 모델링(enzyme-constrained modeling)** 워크플로우입니다. 이전 변경 내역을 포함한 전체 히스토리는 [CHANGELOG.md](CHANGELOG.md)를 참고하세요.

### 🧬 효소 제약 모델링 (GECKO)

- **GECKO 3.0 기반 ecModel 빌드**: Analysis 메뉴 > "Enzyme-Constrained Model (GECKO)…"에서 기존 GEM을 효소 제약 모델(ecModel)로 확장할 수 있습니다. GECKO 3.0 논문 방법론을 충실히 구현하여 Full 및 Light 두 formulation을 모두 지원합니다.

- **통합 워크플로우 다이얼로그**: 사이드바 네비게이션 기반의 단일 다이얼로그에서 GECKO의 모든 단계를 순차적으로 수행할 수 있습니다.
  - **Build ecModel**: Kcat data 및 MW data 로드, sigma / f / Ptot 파라미터 설정, Full/Light formulation 선택
  - **Enzyme Usage Report**: 효소별 사용량 및 제약 분석
  - **Proteomics Integration**: 단백질체 측정 데이터를 통합하여 개별 효소 용량 제약
  - **kcat What-if Analysis**: kcat 값 변경 시 플럭스 영향을 시뮬레이션

- **GECKO YAML 저장/로드 및 프로젝트 상태 복원**: EcStructure 메타데이터를 포함한 GECKO YAML 파일을 저장하고 불러올 수 있습니다. 또한 `.cna` 프로젝트를 저장/재오픈할 때 ecModel 파라미터와 빌드 결과(단백질 pool, 효소 제약, isozyme 매핑 등)가 자동으로 복원되므로 매번 다시 빌드할 필요가 없습니다. Plain GEM YAML과 ecModel YAML을 모두 지원합니다.

- **File → New project from YAML**: GECKO YAML 파일을 직접 불러와 새 프로젝트를 바로 시작할 수 있습니다.

- **번들 GECKO 예제 데이터셋**: File > "New project from GECKO example" 메뉴에서 **Human-GEM**과 **Yeast-GEM** 번들 데이터셋(DLKcat, customKcats, UniProt 포함)으로 GECKO 워크플로우를 즉시 실행할 수 있습니다.

*참고: 이 변경 사항들은 Apache License 2.0 하에 배포되며, 원본 CNApy 프로젝트의 일부입니다.*
