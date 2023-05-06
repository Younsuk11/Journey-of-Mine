## [Machine Learning Project] : Regression

#### Team : 15기 김지호, 염윤석, 우명진, 이제윤
### 영화 관객수 예측
#### link : https://dacon.io/competitions/open/235536/overview/description


### **0.ML library : scikit-learn**

### **1. DATA 구성**
   * title
   * distributor, genre, relase_time, tiem, screening_rat
   * director, dir_prev_bfnum, dir_prev_num, num_staff, num_actor, box_off_num
    
### **2. DATA Preprocessing**
  * Feature selection : 불필요한 변수 삭제 및 자료형 변환
  * 결측치 처리 : 0으로 대체
  * 범주형 변수 : [독점 배급사 5개 + 기타] 로 범주화 | 작품 수 5단위로 묶어서 범주화
  * 개봉 일자 : 성수기/비성수기로 범주화
  
  전처리 사용자 함수를 만들어서 일괄적으로 처리

### **3. Modeling**
  * Lasso Regression
  * Ridge Regression
  * Random Forest Regressor
  * Extra Trees Regressor
  
  ==> GridSearchCV를 통해서 모델 튜닝 후 성능 비교
  
  ==> Voting과 Stacking Ensemble 기법 역시 시도하여 성능 비교

