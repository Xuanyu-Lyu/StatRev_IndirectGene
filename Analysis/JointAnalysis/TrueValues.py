# record the true values of IGE proportions for different simulation conditions
import numpy as np

true_values_pure_ige = {}
true_values_final_ige = {}

f_51 = .3
delta_51 = np.sqrt(0.7885)
v_pure_51 = f_51*delta_51*2
vy_51 = 2.0188
p_pure_51 = v_pure_51 / vy_51

delta_final_51 = np.sqrt(0.7863)
v_final_51 = f_51*delta_final_51*2
vy_final_51 = 2.0188
p_final_51 = v_final_51 / vy_final_51

true_values_pure_ige['05_t1pheVTnoAM'] = p_pure_51
true_values_final_ige['05_t1pheVTnoAM'] = p_final_51

f_52 = 0
delta_52 = np.sqrt(0.8021)
v_pure_52 = f_52*delta_52*2
vy_52 = 1.1977
p_pure_52 = v_pure_52 / vy_52   

delta_final_52 = np.sqrt(0.8001)
v_final_52 = f_52*delta_final_52*2
vy_final_52 = 1.1977
p_final_52 = v_final_52 / vy_final_52

true_values_pure_ige['05_t2socVTnoAM'] = p_pure_52
true_values_final_ige['05_t2socVTnoAM'] = p_final_52

f_61 = 0
delta_61 = np.sqrt(0.8029)
v_pure_61 = f_61*delta_61*2
vy_61 = 1.2664
p_pure_61 = v_pure_61 / vy_61

delta_final_61 = np.sqrt(1.0663)
v_final_61 = f_61*delta_final_61*2
vy_final_61 = 1.2664
p_final_61 = v_final_61 / vy_final_61

true_values_pure_ige['06_t1noVTpheAM'] = p_pure_61
true_values_final_ige['06_t1noVTpheAM'] = p_final_61

f_62 = .3
delta_62 = np.sqrt(0.8029)
v_pure_62 = f_62*delta_62*2
vy_62 = 3.4067
p_pure_62 = v_pure_62 / vy_62
delta_final_62 = np.sqrt(1.0625)
v_final_62 = f_62*delta_final_62*2
vy_final_62 = 3.4067
p_final_62 = v_final_62 / vy_final_62

true_values_pure_ige['06_t2pheVTpheAM'] = p_pure_62
true_values_final_ige['06_t2pheVTpheAM'] = p_final_62

f_71 = 0
delta_71 = np.sqrt(0.7988)
v_pure_71 = f_71*delta_71*2
vy_71 = 0.9990
p_pure_71 = v_pure_71 / vy_71
delta_final_71 = np.sqrt(0.7991)
v_final_71 = f_71*delta_final_71*2
vy_final_71 = 0.9990
p_final_71 = v_final_71 / vy_final_71

true_values_pure_ige['07_t1noVTnoAM'] = p_pure_71
true_values_final_ige['07_t1noVTnoAM'] = p_final_71

f_72 = .3
delta_72 = np.sqrt(0.7993)
v_pure_72 = f_72*delta_72*2
vy_72 = 2.6163
p_pure_72 = v_pure_72 / vy_72
delta_final_72 = np.sqrt(0.8850)
v_final_72 = f_72*delta_final_72*2
vy_final_72 = 2.6163
p_final_72 = v_final_72 / vy_final_72

true_values_pure_ige['07_t2pheVTsocAM'] = p_pure_72
true_values_final_ige['07_t2pheVTsocAM'] = p_final_72

f_81 = 0
delta_81 = np.sqrt(0.8078)
v_pure_81 = f_81*delta_81*2
vy_81 = 1.1756
p_pure_81 = v_pure_81 / vy_81
delta_final_81 = np.sqrt(0.9758)
v_final_81 = f_81*delta_final_81*2
vy_final_81 = 1.1756
p_final_81 = v_final_81 / vy_final_81

true_values_pure_ige['08_t1noVTgenAM'] = p_pure_81
true_values_final_ige['08_t1noVTgenAM'] = p_final_81

f_82 = .3
delta_82 = np.sqrt(0.8106)
v_pure_82 = f_82*delta_82*2
vy_82 = 2.8354
p_pure_82 = v_pure_82 / vy_82
delta_final_82 = np.sqrt(0.9787)
v_final_82 = f_82*delta_final_82*2
vy_final_82 = 2.8354
p_final_82 = v_final_82 / vy_final_82

true_values_pure_ige['08_t2pheVTgenAM'] = p_pure_82
true_values_final_ige['08_t2pheVTgenAM'] = p_final_82

# print the true values
print("True values of pure IGE proportions:")
for condition, value in true_values_pure_ige.items():
    print(f"{condition}: {value:.4f}")
print("\nTrue values of final IGE proportions:")
for condition, value in true_values_final_ige.items():
    print(f"{condition}: {value:.4f}")
