## Deep Learning NLP Project

#### Team : 15기 고태영, 반민정, 염윤석, 최경석
### 병원 개/폐업 분류 예측 경진대회
#### link : https://dacon.io/competitions/official/9565/overview/description

"병원 재무 데이터와 AI를 활용하여, 병원의 개업 또는 폐업 여부를 분류 예측"

### **0. DL library : pytorch**

### **1. DATA 구성**
   * Target : open / close
   * features
      * Sido : 대한민국 행정구역(시,도) - 16곳의 광역 지역
      * Inst Kind : 병원의 종류 - 7가지
      * Owner Change : 대표자 변동여부
      * Bed Count : 병원이 가지고 있는 병상의 수
      * Open Date : 병원의 개원 날짜
      * Income State : revenue | salecost | sga | salary | noi | noe | Interest | Ctax | Profit
      * Asset : liquiedAsset | quickAsset | receivableS | inventoryAsset | nonCAsset | tanAsset | OnonCAsset | receivableL
      * Debt : debt | liquidLiabilities | shortLoan | NCLiabilities | longLoan
      * Equity : netAsset | surplus
    
### **2. DATA Preprocessing**
  * Feature selection : 불필요한 변수 삭제 및 자료형 변환, 개원날짜 => 운영기간 변수 생성
  * 결측치 처리
  * Data Scaling & Log1p Scaling
  * Categorical variable => Dummy Variable
  * Creatig New variables

### **3. Modeling**
  * RandomForestClassifier
  * KNeighborsClassifier
  * ExtraTreeClassifier
  * RidgeClassfier
  * LGBClassfier
  * XGBClassifier
  * GradientBoostigClassifier
  * CatbookstClassifier
  
  ==> GridSearchCV를 통해서 8가지 모델 튜닝 후 성능 비교


### **RESULT : Accuracy**
**ExtraTreeClassififer** 
    
    : public 114 ACC :0.85714
    
**AutoML-Pycaret** : Top5 모델에 대하여 stacking ensemble 적용
top4 : KNN, Dummy classifier, RandomForest, ExtraTree, CatBoost

    : ACC : 0.8888

