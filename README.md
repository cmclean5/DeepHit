# DeepSurv Model
R implementation of python [DeepSurv competing risk model](https://github.com/cmclean5/PublicHealthModels/issues/18)

Paper proposes two step approach to model survival data using deep learning. In the first step, we transform each subject's survival time into a series of jackknife pseudo conditional survival probabilities and then use these pseudo probabilities as a quantitative response variable in the deep neural network model. By using the pseudo values, we reduce a complex survival analysis to a standard regression problem, which greatly simplifies the neural network construction. Our two-step approach is simple, yet very flexible in making risk predictions for survival data, which is very appealing from the practice point of view.

---

