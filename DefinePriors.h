#ifndef PRIOR_H
#define PRIOR_H

#define ADULT_MOR_X0_MEAN 3.8  // according to Keesing 1992, ori:45
#define ADULT_MOR_X0_STD 0.1 // orignal 5
#define ADULT_MOR_K_MEAN -2.30 // according to Keesing 1992, ori:0.1
#define ADULT_MOR_K_STD 0.12 // orignal 0.02
#define BETA_SIG_MEAN 1 // lognormal - normal 0 
#define BETA_SIG_STD 0.2// lognormal - normal 2
#define BETA_MIN 0
#define BETA_MAX 0.4
#define LAR_DAILY_MOR_MEAN -1.93 // Monthly larval mortality, mortality during metamorphism is 85% and daily mortalit is 6%
#define LAR_DAILY_MOR_STD 0.2
#define JUV_SUR_MEAN -2 //Monthly,  after recruitment, lognormal, orignal 0.4
#define JUV_SUR_STD 0.2//Monthly, orignal 0.05
#define DETECT_LO 0.0001
#define DETECT_HI 5
#define RE_MEAN -1.90 
#define RE_STD 0.30
#define DIS_MEAN -3.0
#define DIS_STD 0.3
#define VAR_ERROR_MAX 0.5
#define VAR_ERROR_MIN 0.0001
#endif 
