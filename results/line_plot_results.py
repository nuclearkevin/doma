import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'text.usetex': True})

def error_rel(result, reference):
  return 100.0 * (result - reference) / reference

WIDTH = 16
HEIGHT = 10
PAD = 0.5
LINEWIDTH = 2

kobayashi_1_dims = np.array([5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0])

#-----------------------------------------------------------------------------------------------------------------------------------------------
# 1a: x = z = 5, y \in {5, 15, 25, 35, 45, 55, 65, 75, 85}
#-----------------------------------------------------------------------------------------------------------------------------------------------
kobayashi_1a_bench = np.array([8.29260, 1.87028, 0.713986, 0.384685, 0.253984,  0.137220,  0.0465913,  0.0158766,  0.00547036])

# Diamond differences.
kobayashi_1a_dd_angle_10_space_1 = np.array([7.19185, 1.99889, 0.55205,  0.387174, 0.23297,  0.145836, 0.0436164, 0.01428,   0.00278381])
kobayashi_1a_dd_angle_10_space_2 = np.array([10.2179, 3.37977, 0.92256,  0.457335, 0.291582, 0.191621, 0.0654322, 0.0216539, 0.00718362])
kobayashi_1a_dd_angle_10_space_3 = np.array([8.59579, 1.8066,  0.736503, 0.392815, 0.23564,  0.120626, 0.0394013, 0.0134349, 0.00571944])
kobayashi_1a_dd_angle_10_space_4 = np.array([9.11955, 2.36546, 0.823292, 0.468751, 0.297835, 0.152063, 0.0469145, 0.016721,  0.00655023])

kobayashi_1a_dd_angle_2_space_4  = np.array([9.37995, 2.60431, 0.828742, 0.299079, 0.0748675, 0.0243906, 0.00789227, 0.00335086, 0.00210399])
kobayashi_1a_dd_angle_4_space_4  = np.array([9.13307, 2.4326,  0.83956,  0.34987,  0.303303,  0.170693,  0.0400915,  0.0100397,  0.00309724])
kobayashi_1a_dd_angle_6_space_4  = np.array([9.08546, 2.38337, 0.921323, 0.404173, 0.213837,  0.151771,  0.0595521,  0.0197149,  0.00564349])
kobayashi_1a_dd_angle_8_space_4  = np.array([9.10536, 2.36284, 0.865336, 0.47699,  0.254656,  0.132536,  0.0493567,  0.0192463,  0.00722888])

# Step characteristics.
kobayashi_1a_sc_angle_2_space_4  = np.array([8.5292,  2.34868, 0.868642, 0.393162, 0.181814, 0.0671682, 0.0210858, 0.00770661, 0.00315376])
kobayashi_1a_sc_angle_4_space_4  = np.array([8.5332,  2.30838, 0.930711, 0.490975, 0.303727, 0.166061,  0.0591798, 0.0208698,  0.00743639])
kobayashi_1a_sc_angle_6_space_4  = np.array([8.53417, 2.30374, 0.934177, 0.499821, 0.317129, 0.180716,  0.0671726, 0.024804,   0.00920516])
kobayashi_1a_sc_angle_8_space_4  = np.array([8.53469, 2.30272, 0.934227, 0.501132, 0.319987, 0.183957,  0.0688844, 0.0256666,  0.00963336])
kobayashi_1a_sc_angle_10_space_4 = np.array([8.535,   2.30245, 0.934069, 0.501197, 0.320624, 0.184949,  0.0694544, 0.0259525,  0.00976988])

kobayashi_1a_sc_angle_10_space_1 = np.array([7.45334, 1.78905, 0.983666, 0.614122, 0.420649, 0.222729, 0.101326,  0.0462193, 0.0211585])
kobayashi_1a_sc_angle_10_space_2 = np.array([8.58388, 2.57154, 1.12449,  0.616915, 0.391955, 0.237463, 0.0966251, 0.0390392, 0.0158635])
kobayashi_1a_sc_angle_10_space_3 = np.array([7.75354, 1.88425, 0.845863, 0.477561, 0.313713, 0.167409, 0.0644078, 0.0247834, 0.00959988])

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1a - Diamond Difference Spatial Refinement")
ax1.set_xlabel("$x = 5$, $z = 5$, $y \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_10_space_1, kobayashi_1a_bench),  label='$20^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_10_space_2, kobayashi_1a_bench),  label='$40^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_10_space_3, kobayashi_1a_bench),  label='$60^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_10_space_4, kobayashi_1a_bench),  label='$80^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1a_spatial_refinement_dd.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1a - Step Characteristics Spatial Refinement")
ax1.set_xlabel("$x = 5$, $z = 5$, $y \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_10_space_1, kobayashi_1a_bench),  label='$20^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_10_space_2, kobayashi_1a_bench),  label='$40^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_10_space_3, kobayashi_1a_bench),  label='$60^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_10_space_4, kobayashi_1a_bench),  label='$80^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1a_spatial_refinement_sc.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1a - Diamond Difference Angular Refinement")
ax1.set_xlabel("$x = 5$, $z = 5$, $y \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_2_space_4, kobayashi_1a_bench),  label='$32$ Angles ($n_l, n_c = 2$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_4_space_4, kobayashi_1a_bench),  label='$128$ Angles ($n_l, n_c = 4$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_6_space_4, kobayashi_1a_bench),  label='$288$ Angles ($n_l, n_c = 6$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_8_space_4, kobayashi_1a_bench),  label='$512$ Angles ($n_l, n_c = 8$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_dd_angle_10_space_4, kobayashi_1a_bench), label='$800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.38))
plt.savefig("./line_plots/1a_angular_refinement_dd.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1a - Step Characteristics Angular Refinement")
ax1.set_xlabel("$x = 5$, $z = 5$, $y \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_2_space_4, kobayashi_1a_bench),  label='$32$ Angles ($n_l, n_c = 2$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_4_space_4, kobayashi_1a_bench),  label='$128$ Angles ($n_l, n_c = 4$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_6_space_4, kobayashi_1a_bench),  label='$288$ Angles ($n_l, n_c = 6$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_8_space_4, kobayashi_1a_bench),  label='$512$ Angles ($n_l, n_c = 8$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1a_sc_angle_10_space_4, kobayashi_1a_bench), label='$800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.38))
plt.savefig("./line_plots/1a_angular_refinement_sc.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1a - Fluxes")
ax1.set_xlabel("$x = 5$, $z = 5$, $y \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Scalar Neutron Flux ($cm^{-2}s^{-1}$)")
ax1.plot(kobayashi_1_dims, kobayashi_1a_sc_angle_10_space_4, label='Step Characteristics: $800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, kobayashi_1a_dd_angle_10_space_4, label='Diamond Differences: $800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, kobayashi_1a_bench, label='Reference: Grey Monte Carlo', linewidth=LINEWIDTH)
ax1.set_yscale('log')

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.28))
plt.savefig("./line_plots/1a_results.png", format='png')
plt.show()

#-----------------------------------------------------------------------------------------------------------------------------------------------
# 1b: x, y, z \in {5, 15, 25, 35, 45, 55, 65, 75, 85, 95}
#-----------------------------------------------------------------------------------------------------------------------------------------------
kobayashi_1b_bench = np.array([8.29260, 0.663233, 0.268828, 0.156683, 0.104405, 0.0302145, 0.00406555, 0.000586124, 8.66059e-5])

# Diamond differences.
kobayashi_1b_dd_angle_10_space_1 = np.array([7.19185, 0.715966, 0.264954, 0.160305, 0.0991094, 0.0253011, 0.0052335,  0.000647166, 0.000121009])
kobayashi_1b_dd_angle_10_space_2 = np.array([10.2179, 0.972768, 0.340734, 0.156377, 0.132868,  0.0450985, 0.00703521, 0.000837118, 0.000194821])
kobayashi_1b_dd_angle_10_space_3 = np.array([8.59579, 0.651424, 0.265047, 0.164646, 0.104503,  0.0243088, 0.0036571,  0.000742012, 0.000123662])
kobayashi_1b_dd_angle_10_space_4 = np.array([9.11955, 0.781507, 0.270655, 0.164834, 0.104798,  0.0334191, 0.00574311, 0.000883106, 9.10411e-05])

kobayashi_1b_dd_angle_2_space_4  = np.array([9.37995, 1.12666,  0.105049, 0.0528371, 0.0438104, 0.00489116, 0.00129752, 2.71851e-05,  -2.12835e-05])
kobayashi_1b_dd_angle_4_space_4  = np.array([9.13307, 0.866375, 0.314872, 0.17313,   0.114543,  0.0306458,  0.00176476, -0.000256009, -5.37467e-05])
kobayashi_1b_dd_angle_6_space_4  = np.array([9.08546, 0.798022, 0.327945, 0.145621,  0.116123,  0.0491029,  0.00794364, 0.00080152,    0.0001276])
kobayashi_1b_dd_angle_8_space_4  = np.array([9.10536, 0.780249, 0.282965, 0.169819,  0.118465,  0.0391823,  0.00650347, 0.000859389,   9.89639e-05])

# Step characteristics.
kobayashi_1b_sc_angle_10_space_1 = np.array([7.45334, 0.322224, 0.158877, 0.0984482, 0.0667116, 0.0108549, 0.00210086, 0.000440889, 9.35397e-05])
kobayashi_1b_sc_angle_10_space_2 = np.array([8.58388, 0.599183, 0.237774, 0.137015,  0.0903323, 0.0254387, 0.00394583, 0.000694261, 0.000126714])
kobayashi_1b_sc_angle_10_space_3 = np.array([7.75354, 0.504064, 0.221927, 0.132544,  0.0884854, 0.0182289, 0.00277909, 0.000466385, 7.98498e-05])
kobayashi_1b_sc_angle_10_space_4 = np.array([8.535,   0.625178, 0.252163, 0.145975,  0.0966731, 0.0256637, 0.0037105,  0.000597692, 9.92977e-05])

kobayashi_1b_sc_angle_2_space_4  = np.array([8.5292,  0.800427, 0.223366, 0.0856713, 0.0384689, 0.00583108, 0.00083276, 0.000139178, 2.29644e-05])
kobayashi_1b_sc_angle_4_space_4  = np.array([8.5332,  0.65981,  0.243583, 0.147477,  0.102387,  0.0273763,  0.00392007, 0.000618059, 0.000100138])
kobayashi_1b_sc_angle_6_space_4  = np.array([8.53417, 0.634163, 0.251324, 0.146322,  0.0968938, 0.0255385,  0.00371131, 0.000601796, 0.000100704])
kobayashi_1b_sc_angle_8_space_4  = np.array([8.53469, 0.62745,  0.252127, 0.145954,  0.0967644, 0.0256671,  0.00371226, 0.000597548, 9.91559e-05])

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1b - Diamond Difference Spatial Refinement")
ax1.set_xlabel("$x,y,z \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_10_space_1, kobayashi_1b_bench),  label='$20^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_10_space_2, kobayashi_1b_bench),  label='$40^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_10_space_3, kobayashi_1b_bench),  label='$60^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_10_space_4, kobayashi_1b_bench),  label='$80^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1b_spatial_refinement_dd.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1b - Step Characteristic Spatial Refinement")
ax1.set_xlabel("$x,y,z \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_10_space_1, kobayashi_1b_bench),  label='$20^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_10_space_2, kobayashi_1b_bench),  label='$40^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_10_space_3, kobayashi_1b_bench),  label='$60^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_10_space_4, kobayashi_1b_bench),  label='$80^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1b_spatial_refinement_sc.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1b - Diamond Difference Angular Refinement")
ax1.set_xlabel("$x,y,z \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_2_space_4, kobayashi_1b_bench),  label='$32$ Angles ($n_l, n_c = 2$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_4_space_4, kobayashi_1b_bench),  label='$128$ Angles ($n_l, n_c = 4$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_6_space_4, kobayashi_1b_bench),  label='$288$ Angles ($n_l, n_c = 6$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_8_space_4, kobayashi_1b_bench),  label='$512$ Angles ($n_l, n_c = 8$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_dd_angle_10_space_4, kobayashi_1b_bench), label='$800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.38))
plt.savefig("./line_plots/1b_angular_refinement_dd.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1b - Step Characteristics Angular Refinement")
ax1.set_xlabel("$x,y,z \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_2_space_4, kobayashi_1b_bench),  label='$32$ Angles ($n_l, n_c = 2$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_4_space_4, kobayashi_1b_bench),  label='$128$ Angles ($n_l, n_c = 4$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_6_space_4, kobayashi_1b_bench),  label='$288$ Angles ($n_l, n_c = 6$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_8_space_4, kobayashi_1b_bench),  label='$512$ Angles ($n_l, n_c = 8$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1b_sc_angle_10_space_4, kobayashi_1b_bench), label='$800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.38))
plt.savefig("./line_plots/1b_angular_refinement_sc.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1b - Fluxes")
ax1.set_xlabel("$x,y,z \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Scalar Neutron Flux ($cm^{-2}s^{-1}$)")
ax1.plot(kobayashi_1_dims, kobayashi_1b_sc_angle_10_space_4, label='Step Characteristics: $800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, kobayashi_1b_dd_angle_10_space_4, label='Diamond Differences: $800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, kobayashi_1b_bench, label='Reference: Grey Monte Carlo', linewidth=LINEWIDTH)
ax1.set_yscale('log')

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.28))
plt.savefig("./line_plots/1b_results.png", format='png')
plt.show()

#-----------------------------------------------------------------------------------------------------------------------------------------------
# 1c: y = 55, z = 5, x \in {5, 15, 25, 35, 45, 55, 65, 75, 85, 95}
#-----------------------------------------------------------------------------------------------------------------------------------------------
kobayashi_1c_bench = np.array([0.137220, 0.127890, 0.113582, 0.0959578, 0.0782701, 0.0567030, 0.0188631,  0.00646624,  0.00228099])

# Diamond differences.
kobayashi_1c_dd_angle_10_space_1 = np.array([0.145836, 0.137323, 0.118353, 0.107105,  0.0868175, 0.0562456, 0.0217341, 0.00327227, 0.00310898])
kobayashi_1c_dd_angle_10_space_2 = np.array([0.191621, 0.184397, 0.148391, 0.135251,  0.130069,  0.078181,  0.0222391, 0.00935347, 0.00342325])
kobayashi_1c_dd_angle_10_space_3 = np.array([0.120626, 0.117482, 0.110007, 0.0957292, 0.073202,  0.046566,  0.0161281, 0.00743358, 0.00281286])
kobayashi_1c_dd_angle_10_space_4 = np.array([0.152063, 0.150671, 0.137018, 0.118466,  0.101546,  0.0596483, 0.0227721, 0.00828055, 0.00228161])

kobayashi_1c_dd_angle_2_space_4  = np.array([0.0243906, 0.0241489, 0.0359937, 0.0324809, 0.0172694, 0.0115145, 0.00593909, 0.00618112, 0.00211057])
kobayashi_1c_dd_angle_4_space_4  = np.array([0.170693,  0.161568,  0.117316,  0.131916,  0.0778481, 0.0473661, 0.0322303,  0.00975141, 0.00245024])
kobayashi_1c_dd_angle_6_space_4  = np.array([0.151771,  0.16727,   0.104379,  0.121062,  0.0940065, 0.0690771, 0.0237085,  0.00741201, 0.0018982])
kobayashi_1c_dd_angle_8_space_4  = np.array([0.132536,  0.143906,  0.133644,  0.108427,  0.0812053, 0.0676862, 0.0184413,  0.00795957, 0.00241719])

# Step characteristics.
kobayashi_1c_sc_angle_10_space_1 = np.array([0.222729, 0.134532, 0.109248, 0.0890307, 0.0719821, 0.0414663, 0.0211041, 0.0104485,  0.00510654])
kobayashi_1c_sc_angle_10_space_2 = np.array([0.237463, 0.179115, 0.148792, 0.123735,  0.101315,  0.0669976, 0.0296846, 0.0125914,  0.00533747])
kobayashi_1c_sc_angle_10_space_3 = np.array([0.167409, 0.133901, 0.113071, 0.0939785, 0.0762653, 0.0476668, 0.0205365, 0.0083451,  0.00337837])
kobayashi_1c_sc_angle_10_space_4 = np.array([0.184949, 0.154974, 0.1318,   0.110468,  0.0904549, 0.0599194, 0.0244681, 0.00944669, 0.00368219])

kobayashi_1c_sc_angle_2_space_4  = np.array([0.0671682, 0.101917, 0.0858999, 0.0445408, 0.0288913, 0.0243044, 0.0154018, 0.00783347, 0.00354548])
kobayashi_1c_sc_angle_4_space_4  = np.array([0.166061,  0.1591,   0.119177,  0.102776,  0.0854703, 0.0554903, 0.0236347, 0.00932802, 0.00354729])
kobayashi_1c_sc_angle_6_space_4  = np.array([0.180716,  0.159009, 0.128275,  0.10928,   0.0892208, 0.0590476, 0.024146,  0.00933467, 0.0036496])
kobayashi_1c_sc_angle_8_space_4  = np.array([0.183957,  0.156573, 0.131025,  0.110136,  0.090179,  0.0597112, 0.0243967, 0.00942225, 0.00367276])

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1c - Diamond Difference Spatial Refinement")
ax1.set_xlabel("$y = 55$, $z = 5$, $x \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_10_space_1, kobayashi_1c_bench),  label='$20^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_10_space_2, kobayashi_1c_bench),  label='$40^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_10_space_3, kobayashi_1c_bench),  label='$60^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_10_space_4, kobayashi_1c_bench),  label='$80^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.33))
plt.savefig("./line_plots/1c_spatial_refinement_dd.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1c - Step Characteristic Spatial Refinement")
ax1.set_xlabel("$y = 55$, $z = 5$, $x \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_10_space_1, kobayashi_1c_bench),  label='$20^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_10_space_2, kobayashi_1c_bench),  label='$40^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_10_space_3, kobayashi_1c_bench),  label='$60^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_10_space_4, kobayashi_1c_bench),  label='$80^3$ Cells, $800$ Angles', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1c_spatial_refinement_sc.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1c - Diamond Difference Angular Refinement")
ax1.set_xlabel("$y = 55$, $z = 5$, $x \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_2_space_4, kobayashi_1c_bench),  label='$32$ Angles ($n_l, n_c = 2$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_4_space_4, kobayashi_1c_bench),  label='$128$ Angles ($n_l, n_c = 4$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_6_space_4, kobayashi_1c_bench),  label='$288$ Angles ($n_l, n_c = 6$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_8_space_4, kobayashi_1c_bench),  label='$512$ Angles ($n_l, n_c = 8$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_dd_angle_10_space_4, kobayashi_1c_bench), label='$800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1c_angular_refinement_dd.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1c - Step Characteristics Angular Refinement")
ax1.set_xlabel("$y = 55$, $z = 5$, $x \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Relative Error ($\%$)")
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_2_space_4, kobayashi_1c_bench),  label='$32$ Angles ($n_l, n_c = 2$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_4_space_4, kobayashi_1c_bench),  label='$128$ Angles ($n_l, n_c = 4$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_6_space_4, kobayashi_1c_bench),  label='$288$ Angles ($n_l, n_c = 6$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_8_space_4, kobayashi_1c_bench),  label='$512$ Angles ($n_l, n_c = 8$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, error_rel(kobayashi_1c_sc_angle_10_space_4, kobayashi_1c_bench), label='$800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.87))
plt.savefig("./line_plots/1c_angular_refinement_sc.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("3D Kobayashi Benchmark 1c - Fluxes")
ax1.set_xlabel("$y = 55$, $z = 5$, $x \in [5, 85]$ ($cm$)")
ax1.set_ylabel("Scalar Neutron Flux ($cm^{-2}s^{-1}$)")
ax1.plot(kobayashi_1_dims, kobayashi_1c_sc_angle_10_space_4, label='Step Characteristics: $800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, kobayashi_1c_dd_angle_10_space_4, label='Diamond Differences: $800$ Angles ($n_l, n_c = 10$), $80^3$ Cells', linewidth=LINEWIDTH)
ax1.plot(kobayashi_1_dims, kobayashi_1c_bench, label='Reference: Grey Monte Carlo', linewidth=LINEWIDTH)
ax1.set_yscale('log')

fig.legend(loc='upper left', bbox_to_anchor=(0.13, 0.28))
plt.savefig("./line_plots/1c_results.png", format='png')
plt.show()
