# Head-Neck-Cancer--ROC--SVM--KM--Expression-Analysis collaboration with Dr. Stanam, University of Iowa, Iowa City, USA. https://www.linkedin.com/in/aditya-stanam-07bbb819/
Springer Journal of Maxillofacial and Oral Surgery. [Paper Link:](https://www.ncbi.nlm.nih.gov/pubmed/30996559)
[R: survival, Biobase, pbapply]. Date created/updated: December, 9, 2024.
<img width="749" alt="Six gene" src="https://github.com/spawar2/Head-Neck-Cancer--ROC--SVM--KM--Expression-Analysis/assets/25118302/085134c5-2fa0-48d0-acc0-f13a33600862">

<img width="738" alt="Screenshot 2023-05-17 at 10 26 55 AM" src="https://github.com/spawar2/Head-Neck-Cancer--ROC--SVM--KM--Expression-Analysis/assets/25118302/e84512a9-b76b-446a-bcab-30f06d36fd61">

Georgia State University, Atlanta, United States of America (USA).
https://csds.gsu.edu/

ROC.R: Head-Neck-Cancer microarray data read, KM COX analysis was performed on survial and days, with significant P value, Survival Kaplan Meir KM analysis, area under curve and receiver operating characteristic, AUC-ROC evaluation.
selected function(re.search,getROC_AUC, plot, as.numeric, sapply, as.numeric, unlist, svm, predict, Surv, survfit, surv_pvalue).
A Six-Gene-Based Prognostic Model Predicts Survival in Head and Neck Squamous Cell Carcinoma Patients, Shrikant Pawar and Stanam, A., Publication: Springer: Journal of Maxillofacial and Oral Surgery (Publication date: January 4), collaboration with Dr. Stanam, University of Iowa, Iowa City, USA, IF=1.0.Github, Article link, [Cited times: 3]^^^ DOI: 10.1007/s12663-019-01187-z, Issue: 2, Volume: 18, Pages: 320-327.
†
†Corresponding author. ††First author. †††Second author. ††††Third author. †††††author.
Testing: Prediction_train
        alive dead
  alive    34    0
  dead      1    1
((34+1)/(nrow(training)))*100
[1] 97.22222.

table(testing$V2,pred_test)
Prediction_test
        alive dead
  alive   214    5
  dead     31   11
((214+11)/(nrow(testing)))*100
[1] 86.2069.
