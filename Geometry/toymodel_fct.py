import matplotlib.pyplot as plt
import math
import numpy as np

class Toy_FCT:
	def __init__(self):
		# Default paramters. Can be changed by calling functions
		# FCT Parameters
		self.z_fct = [4.42, 4.44, 4.46, 4.48, 4.50, 4.52, 4.60, 4.70, 4.80, 4.90, 5.0] # Distance discs to origin (m)
		self.r_out_fct = [0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.18, 0.18, 0.19, 0.19] # Outer radii discs (m)
		self.r_in_fct = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05] # Inner radii discs (m)
		self.layerx2X0_def_fct = 0.01 # Thickness of the layers default value in % of a radiation length. Will fill the list if no values are given
		self.layerx2X0_list_fct = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01] # Thickness of the layers list
		self.layerType_def_fct = 0 # Default type of the layers used for the FCT. 0: Si active disc, 1: Si active square, 2: Pb passive converter disc (see O2 code). Will fill the list if empty
		self.layerType_list_fct = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # Type of the layers list

		# FCT Magnet
		self.z_fct_magnet = 0 # (m)
		self.length_fct_magnet = 0 # (m)
		self.r_in_fct_magnet = 0 # (m)
		self.r_out_fct_magnet = 0 # (m)
		self.set_fct_magnet_values()

		# Iris Parameters
		self.z_its = [0.26, 0.3, 0.34] # (m)
		self.r_out_its = [0.025, 0.025, 0.025] # (m)
		self.r_in_its = [0.005, 0.005, 0.005] # (m)
		self.layerx2X0_def_its = 0.001 # Thickness of the layers default value in % of a radiation length. Will fill the list if no values are given
		self.layerx2X0_list_its = [] # Thickness of the layers list

		# FT3 Parameters - Symmetric disc setup. Will make discs in both the positive and negative direction
		self.z_ft3 = [0.77, 1., 1.22, 1.50, 1.80, 2.20, 2.60, 3.00, 3.50] # Absolute distance to origin (m)
		self.r_out_ft3 = [0.35, 0.35, 0.35, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68] # Outer radii discs (m)
		self.r_in_ft3_A_side = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05] # Inner radii discs (m) - A-side
		self.r_in_ft3_C_side = [0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07] # Inner radii discs (m) - C-side
		self.layerx2X0_def_ft3 = 0.01 # Thickness of the layers default value in % of a radiation length. Will fill the list if no values are given
		self.layerx2X0_list_ft3 = [] # Thickness of the layers list

		# TRK Parameters
		self.z_trk = [0.25, 0.25, 0.25, 0.62, 0.62, 0.62, 0.62, 0.62, 1.29, 1.29, 1.29] # Half-length of the barrel layers (m)
		self.r_trk = [0.005, 0.012, 0.025, 0.07, 0.09, 0.12, 0.2, 0.3, 0.45, 0.6, 0.8] # Radius of the barrel laters (m)

		# Coldplate
		self.z_coldplate = 0.25 # Half-length (m)
		self.r_coldplate = 0.026 # (m)

		# Beam pipe and vacuum chamber Parameters
		self.r_bp_A_side = 0.018 # Radius beampipe A-side (m)
		self.r_bp_C_side = 0.056 # Radius beampipe C-side (m)
		self.z_vacuumVesselWall = 0.42 # (m). Value is a guess based on the schematic of the ALICE 3 Tracker in the LOI, Figure 77
		self.r_vacuumVessel = 0.056 # Inner radius beam pipe its. 4th barrel layer is at r = 0.0375, so this is a guess
		self.eta_max = 5.1
		self.eta_min = 3.8
		self.window = False
		self.window_design = 3
		self.z_windowLength = 1.5 # (m)

		# Services vacuum vessel
		self.z_midpoint_services_vacv = -0.2 # (m). Offset for the vacuum vessel services. Makes it easier to draw
		self.z_services_vacv = self.z_vacuumVesselWall - 0.02 - self.z_midpoint_services_vacv # (m) Half length
		self.r_in_services_vacv = self.r_out_its[2] + 0.003
		self.r_out_services_vacv = self.r_in_services_vacv + 0.01

		# TOF
		self.r_tof = 0.85 # (m)
		self.z_tof = 3.5 # Half-length (m)

		# iTOF
		self.r_itof = 0.19 # (m)
		self.z_itof = 0.62 # Half-length (m)

		# FTOF
		self.z_ftof = 3.7 # (m)
		self.r_in_ftof = 0.15 # (m)
		self.r_out_ftof = 1. # (m)

		# RICH
		self.z_rich = 3.5 # Half-length (m)
		self.r_in_rich = 0.9 # (m)
		self.r_out_rich = 1.2 # (m)

		# Forward RICH
		self.z_frich = 3.75 # (m)
		self.length_frich = 0.4 # (m) Length of fRICH in z direction
		self.r_in_frich = 0.15 # (m)
		self.r_out_frich = 1.15 # (m)

		# FOCAL - Since the FOCAL is not azimuthally symmetric, and is in fact more cube shaped, this just outliens the position where it would roughly be.
		self.z_focal = 7 # (m)
		self.length_focal = 1.5 # (m)
		self.eta_min_focal = 3.2 # (m)
		self.eta_max_focal = 5.8 # (m)
		self.r_in_lowz_focal = self.z_focal * math.tan(self.pseurap_to_angle(self.eta_max_focal))
		self.r_out_lowz_focal = self.z_focal * math.tan(self.pseurap_to_angle(self.eta_min_focal))
		self.r_in_highz_focal = (self.z_focal + self.length_focal) * math.tan(self.pseurap_to_angle(self.eta_max_focal))
		self.r_out_highz_focal = (self.z_focal + self.length_focal) * math.tan(self.pseurap_to_angle(self.eta_min_focal))

		# FCT RICH parameters
		self.fct_rich_radiator_length = 1.6 # (m) Like the dRICH of the ePIC at the EIC
		self.fct_rich_service_length = 0.2 # (m) Estimated space in z required for the services of the FCT RICH
		self.z_fct_rich = self.z_fct[-1] + 0.1 # 0.1 to have the FCT RICH start at a realistic distance from the last layer of the FCT
		self.r_in_fct_rich = 0.05 # (m) 
		self.r_out_fct_rich_lowz = self.z_fct_rich * math.tan(self.pseurap_to_angle(4)) # (m)
		self.r_out_fct_rich_highz = (self.z_fct_rich + self.fct_rich_radiator_length) * math.tan(self.pseurap_to_angle(4)) # (m)

		# Iris Tracker in the open position. This is added to the 
		self.r_open_closed_difference = 0.016 - 0.0048 # (m)

		# Pseudorapidity lines
		# self.pseudorapidity_list = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
		self.pseudorapidity_list = [4, 5]
		self.pseurap_origins = [0.]
		
		# Draw parameters
		self.fig, self.ax = plt.subplots() # Output plot
		self.xleft = -8.6 # Left part of field of vision
		self.xright = 8.6 # Field of vision
		self.yup = 1. # Top part of field ov vision
		self.color_areas = True
		
		self.initialize_values()


	def initialize_values(self):
		if(not self.layerx2X0_list_fct):
			self.layerx2X0_list_fct = [self.layerx2X0_def_fct for a in self.z_fct]
		if(not self.layerType_list_fct):
			self.layerType_list_fct = [self.layerType_def_fct for a in self.z_fct]
		if(not self.layerx2X0_list_ft3):
			self.layerx2X0_list_ft3 = [self.layerx2X0_def_ft3 for a in self.z_ft3]
		if(not self.layerx2X0_list_its):
			self.layerx2X0_list_its = [self.layerx2X0_def_its for a in self.z_its]

	# Converter functions
	def ang_to_pseurap(self, angle): # Angle in radians
		pseurap = -math.log(math.tan(angle/2.))
		return pseurap
		
	def pseurap_to_angle(self, eta):
		theta = 2 * math.atan(math.exp(-eta))
		if(eta > 0):
			return theta
		else:
			return math.pi - theta
	
	# Set radii of FCT and FT3 discs to only extend up to a certain max pseudorapidity
	def set_inner_disc_radii(self, eta_min, which = 0):
		# 0: FCT, 1: FT3 (Excluding Iris)
		if(which == 0):
			r_list = []
			for z in self.z_fct:
				r = z * math.tan(self.pseurap_to_angle(eta_min))
				if(r < 0.05):
					r_list.append(0.05)
				else:
					r_list.append(r)
			self.r_in_fct = r_list
		elif(which == 1):
			r_list = []
			for z in self.z_ft3:
				r = z * math.tan(self.pseurap_to_angle(eta_min))
				if(r < 0.05):
					r_list.append(0.05)
				else:
					r_list.append(r)
			self.r_in_ft3_A_side = r_list
			self.r_in_ft3_C_side = r_list

	def set_outer_disc_radii(self, eta_max, which = 0):
		# 0: FCT, 1: FT3 (Excluding Iris)
		if(which == 0):
			self.r_out_fct = [0.01 * math.ceil(100. * z * math.tan(self.pseurap_to_angle(eta_max))) for z in self.z_fct]
		elif(which == 1):
			self.r_out_ft3 = [z * math.tan(self.pseurap_to_angle(eta_max)) for z in self.z_ft3]
	
	def set_fct_by_first_layer(self, zFirstLayer):
		spacing = [0, 0.02, 0.02, 0.02, 0.02, 0.02, 0.08, 0.10, 0.10, 0.10, 0.10]
		self.z_fct = [zFirstLayer + sum(spacing[:index]) for index, _ in enumerate(spacing)]
		print(self.z_fct)
		self.set_outer_disc_radii(4, 0)
	
	def set_fct_magnet_values(self):
		# Entirely paramatrised based on the values given for the FCT
		self.z_fct_magnet = self.z_fct[0] - 0.20 # (m)
		self.length_fct_magnet = self.z_fct[-1] - self.z_fct[0] + 0.4 # (m)
		self.r_in_fct_magnet = max(self.r_out_fct) + 0.10 # (m)
		self.r_out_fct_magnet = self.r_in_fct_magnet + 0.3 # (m)

	def fct_rich_in_front_of_fct(self, distance_from_ft3=0.1):
		self.z_fct_rich = self.z_ft3[-1] + distance_from_ft3 # 0.1 to have the FCT RICH start at a realistic distance from the last layer of the FT3
		self.r_out_fct_rich_lowz = self.z_fct_rich * math.tan(self.pseurap_to_angle(4)) # (m)
		self.r_out_fct_rich_highz = (self.z_fct_rich + self.fct_rich_radiator_length) * math.tan(self.pseurap_to_angle(4)) # (m)
		self.set_fct_by_first_layer(self.z_fct_rich + self.fct_rich_radiator_length + self.fct_rich_service_length)
	
	def open_iris(self):
		# ITS discs
		self.r_in_its = [self.r_open_closed_difference + r_in for r_in in self.r_in_its]
		self.r_out_its = [self.r_open_closed_difference + r_out for r_out in self.r_out_its]

		# Barrel layers 
		for i in range(0, 3):
			self.r_trk[i] += self.r_open_closed_difference

		# Services, coldplate, vacuum vessel
		self.r_in_services_vacv = self.r_out_its[2] + 0.003
		self.r_out_services_vacv = self.r_in_services_vacv + 0.01

		self.r_coldplate += self.r_open_closed_difference
	
	def print_disc_details(self, which):
		# 0: FCT, 1: FT3 (including Iris), 2:TRK
		# This function does not work as well anymore for FT3 discs since the A and C side are not symmetric anymore.
		# Should be easy to fix

		if(which == 0):
			print("---------------------------------------------")
			print("Displaying details for the FCT")
			print("")
			z = self.z_fct
			r_out = self.r_out_fct
			r_in = self.r_in_fct
		elif(which == 1):
			print("---------------------------------------------")
			print("Displaying details for the FT3")
			print("")
			z = self.z_its + self.z_ft3
			r_out = self.r_out_its + self.r_out_ft3
			r_in = self.r_in_its + self.r_in_ft3_A_side
		
		for i in range(0, len(z)):
			eta_minmax = self.pseurap_range(abs(z[i]), r_out[i], r_in[i])
			print("Layer ", i)
			print("Distance: ", z[i])
			print("Inner radius of disc: ", "%.3f" % r_in[i], "m")
			print("Max pseudorapidity: ","%.3f" % eta_minmax[1], "m")
			print("Outer radius of disc: ", "%.3f" % r_out[i], "m")
			print("Min pseudorapidity: ", "%.3f" % eta_minmax[0])
			print("Area: ", "%.3f" % (math.pi * (r_out[i]*r_out[i] - r_in[i]*r_in[i])), "m^2")
			print("")
			print("------")
			print("")
		print("The collective area of the (disc) layers: ")
		print("%.2f" % self.discs_area(r_out, r_in), "m2")
		print("---------------------------------------------")
		return

	def pseurap_range(self, d, r_out, r_in):
		# Calculates the pseudorapidity range of a disc
		pseurap_min = self.ang_to_pseurap(math.atan(r_out/d))
		pseurap_max = self.ang_to_pseurap(math.atan(r_in/d))
		return pseurap_min, pseurap_max

	def discs_area(self, r_out, r_in):
		area = 0.

		for i in range(len(r_out)):
			area += math.pi * (r_out[i] ** 2. - r_in[i] ** 2.)
		
		return area
	
	def write_layout(self, which, output_name):
		# Creates an output file for the FCT or FT3 such that it can be read in by O2
		# which: 0: FCT, 1: FT3		
		output = ""
		if(output_name[-4:] != '.cfg'):
			output = output_name + '.cfg'
		else:
			output = output_name
		
		if(which == 0):
			print("Writing an output file for the configuration parameters of the FCT")
			with open(output, 'w') as f:
				f.write('#layerType  z_layer  r_in  r_out  Layerx2X0\n')
				for i in range(0, len(self.z_fct)):
					f.write("%.f" % self.layerType_list_fct[i])
					f.write('  ')
					f.write("%.1f" % (100 * self.z_fct[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_in_fct[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_out_fct[i]))
					f.write('  ')
					f.write(str(self.layerx2X0_list_fct[i]))
					f.write('\n')
		if(which == 1):
			print("Writing an output file for the configuration parameters of the FT3")
			with open(output, 'w') as f:
				f.write('#z_layer  r_in  r_out  Layerx2X0\n')
				# Negative layers first - C-side
				for i in range(0, len(self.z_its)):
					f.write("%.1f" % (-100 * self.z_its[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_in_its[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_out_its[i]))
					f.write('  ')
					f.write(str(self.layerx2X0_list_its[i]))
					f.write('\n')					
				for i in range(0, len(self.z_ft3)):
					f.write("%.1f" % (-100 * self.z_ft3[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_in_ft3_C_side[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_out_ft3[i]))
					f.write('  ')
					f.write(str(self.layerx2X0_list_ft3[i]))
					f.write('\n')
				# Positive layers - A-side
				for i in range(0, len(self.z_its)):
					f.write("%.1f" % (100 * self.z_its[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_in_its[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_out_its[i]))
					f.write('  ')
					f.write(str(self.layerx2X0_list_its[i]))
					f.write('\n')
				for i in range(0, len(self.z_ft3)):
					f.write("%.1f" % (100 * self.z_ft3[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_in_ft3_C_side[i]))
					f.write('  ')
					f.write("%.1f" % (100 * self.r_out_ft3[i]))
					f.write('  ')
					f.write(str(self.layerx2X0_list_ft3[i]))
					f.write('\n')

	def show(self):
		self.ax.set_xlabel("Distance z (m)", fontsize=14)
		self.ax.set_ylabel("Radius (m)", fontsize=14)
		self.ax.tick_params(axis='both', labelsize = 14)
		self.ax.set_ylim(0, self.yup)
		self.ax.set_xlim(self.xleft, self.xright)
		self.ax.legend(prop={'size': 14}, ncol=2)
		self.ax.invert_xaxis()
		plt.show()
		self.fig.show()
		return

	def grid(self):
		self.ax.grid(True)
		return

	def draw_all(self):
		self.draw_FCT()
		self.draw_FCT_MAGNET()
		self.draw_FT3()
		self.draw_TRK()
		self.draw_Iris()
		self.draw_RICH()
		self.draw_FRICH()
		self.draw_TOF()
		self.draw_iTOF()
		self.draw_FTOF()
		self.draw_beam_pipe()
		self.draw_FCT_RICH()
		self.draw_FOCAL()
		self.draw_passive_material()
		self.draw_coldplate()
		self.draw_pseu_lines()
		

	def draw_FCT(self):
		for i in range(0, len(self.z_fct)):
			x = [self.z_fct[i], self.z_fct[i]]
			y = [self.r_in_fct[i], self.r_out_fct[i]]
			if(i == 0):
				self.ax.plot(x, y, 'm-', label="FCT")
			else:
				self.ax.plot(x, y, 'm-')
	
	def draw_FCT_MAGNET(self):
		self.z_fct_magnet = self.z_fct[0] - 0.20 # (m)
		self.length_fct_magnet = self.z_fct[-1] - self.z_fct[0] + 0.4 # (m)
		self.r_in_fct_magnet = max(self.r_out_fct) + 0.10 # (m)
		self.r_out_fct_magnet = self.r_in_fct_magnet + 0.3 # (m)
		self.ax.plot([self.z_fct_magnet, self.z_fct_magnet], [self.r_in_fct_magnet, self.r_out_fct_magnet], color='deeppink', label='Dipole')
		self.ax.plot([self.z_fct_magnet, self.z_fct_magnet + self.length_fct_magnet], [self.r_in_fct_magnet, self.r_in_fct_magnet], color='deeppink')
		self.ax.plot([self.z_fct_magnet, self.z_fct_magnet + self.length_fct_magnet], [self.r_out_fct_magnet, self.r_out_fct_magnet], color='deeppink')
		self.ax.plot([self.z_fct_magnet + self.length_fct_magnet, self.z_fct_magnet + self.length_fct_magnet], [self.r_in_fct_magnet, self.r_out_fct_magnet], color='deeppink')

		if(self.color_areas):
			self.ax.fill_between([self.z_fct_magnet, self.z_fct_magnet + self.length_fct_magnet], [self.r_in_fct_magnet, self.r_in_fct_magnet], [self.r_out_fct_magnet, self.r_out_fct_magnet], facecolor='deeppink')
	
	def draw_Iris(self):
		for i in range(0, len(self.z_its)):
			x1 = [self.z_its[i], self.z_its[i]]
			x2 = [-self.z_its[i], -self.z_its[i]]
			y = [self.r_in_its[i], self.r_out_its[i]]
			if(i == 0):
				self.ax.plot(x1, y, 'y-', label="Iris discs")
				self.ax.plot(x2, y, 'y-')
			else:
				self.ax.plot(x1, y, 'y-')
				self.ax.plot(x2, y, 'y-')	

	def draw_FT3(self):
		# Draw A-side (pos z) discs
		for i in range(0, len(self.z_ft3)):
			x = [self.z_ft3[i], self.z_ft3[i]]
			y = [self.r_in_ft3_A_side[i], self.r_out_ft3[i]]
			if(i == 0):
				self.ax.plot(x, y, 'g-', label="Tracker discs")
			else:
				self.ax.plot(x, y, 'g-')
		
		# Draw C-side (neg z) discs
		for i in range(0, len(self.z_ft3)):
			x = [-self.z_ft3[i], -self.z_ft3[i]]
			y = [self.r_in_ft3_C_side[i], self.r_out_ft3[i]]
			self.ax.plot(x, y, 'g-')
	
	def draw_TRK(self):
		for i in range(0, len(self.r_trk)):
			x = [-self.z_trk[i], self.z_trk[i]]
			y = [self.r_trk[i], self.r_trk[i]]
			if(i == 0):
				self.ax.plot(x, y, 'c-', linewidth=2, label="Barrel layers")
			else:
				self.ax.plot(x, y, 'c-', linewidth=2)
	
	def draw_beam_pipe(self):
		# plot beam line
		self.ax.plot([-self.xright, self.xright], [0., 0.], 'k-')

		# plot beam pipe - C-side
		self.ax.plot([-self.xright, -self.z_vacuumVesselWall], [self.r_bp_C_side, self.r_bp_C_side], 'b-')
		# vacuum vessel - top part
		self.ax.plot([-self.z_vacuumVesselWall, self.z_vacuumVesselWall], [self.r_vacuumVessel, self.r_vacuumVessel], 'b-', label = "Beam pipe")

		if(self.window):
			if(self.window_design == 1):
				# Goes down, then makes a double cone shape. Not favored
				theta_eta_max = self.pseurap_to_angle(self.eta_max)
				theta_eta_min = self.pseurap_to_angle(self.eta_min)
				z1 = self.r_bp_A_side / (math.tan(theta_eta_max))
				z2 = self.r_in_ft3_A_side[0] / (math.tan(theta_eta_min))
				z3 = self.r_bp_A_side / (math.tan(theta_eta_min))
				# Plot pipe 1
				self.ax.plot([self.z_vacuumVesselWall, z3], [self.r_bp_A_side, self.r_bp_A_side], 'b-')
				# Plot Window 1
				self.ax.plot([z2, z1], [self.r_in_ft3_A_side[0], self.r_bp_A_side], 'b-')
				# Plot Window 2
				self.ax.plot([z3, z2], [self.r_bp_A_side, self.r_in_ft3_A_side[0]], 'b-')
				# Plot pipe 2
				self.ax.plot([self.xright, z1], [self.r_bp_A_side, self.r_bp_A_side], 'b-')
				# Plot vacuum vessel wall - pos z side
				self.ax.plot([self.z_vacuumVesselWall, self.z_vacuumVesselWall], [self.r_bp_A_side, self.r_vacuumVessel], 'b-')
			elif(self.window_design == 2):
				# Elongates the primary vacuum vessel and then makes a wall. Intersects iwth the FT3 discs
				theta_eta_max = self.pseurap_to_angle(self.eta_max)
				theta_eta_min = self.pseurap_to_angle(self.eta_min)

				vacVExt = self.r_vacuumVessel / math.tan(theta_eta_min) # Absolute value
				win = self.r_bp_A_side / math.tan(theta_eta_max) # Absolute value

				# Plot vacuum vessel extension
				self.ax.plot([self.z_vacuumVesselWall, vacVExt], [self.r_vacuumVessel, self.r_vacuumVessel], 'b-')
				# Plot Window
				self.ax.plot([win, vacVExt], [self.r_bp_A_side, self.r_vacuumVessel], 'b-')
				# Plot pipe A side
				self.ax.plot([win, self.xright], [self.r_bp_A_side, self.r_bp_A_side], 'b-')
			elif(self.window_design == 3):
				# From the primary vacuum vessel, makes a cone that then connects down to the beam pipe
				# Plot window
				self.ax.plot([self.z_vacuumVesselWall +  self.z_windowLength, self.z_vacuumVesselWall], [self.r_bp_A_side , self.r_bp_C_side], 'b-')
				# Plot pipe A side
				self.ax.plot([self.xright, self.z_vacuumVesselWall + self.z_windowLength], [self.r_bp_A_side, self.r_bp_A_side], 'b-')

		else:
			self.ax.plot([self.z_vacuumVesselWall, self.z_vacuumVesselWall], [self.r_bp_A_side, self.r_vacuumVessel], 'b-')
			self.ax.plot([self.z_vacuumVesselWall, self.xright], [self.r_bp_A_side, self.r_bp_A_side], 'b-')
	
	def draw_TOF(self):
		self.ax.plot([-self.z_tof, self.z_tof], [self.r_tof, self.r_tof], color='peru', label='TOF')
	
	def draw_iTOF(self):
		self.ax.plot([-self.z_itof, self.z_itof], [self.r_itof, self.r_itof], color='saddlebrown', label='iTOF')

	def draw_FTOF(self):
		self.ax.plot([self.z_ftof, self.z_ftof], [self.r_in_ftof, self.r_out_ftof], color='tab:orange', label='fTOF')
		self.ax.plot([-self.z_ftof, -self.z_ftof], [self.r_in_ftof, self.r_out_ftof], color='tab:orange')

	def draw_RICH(self):
		self.ax.plot([-self.z_rich, self.z_rich], [self.r_in_rich, self.r_in_rich], color='lightsteelblue', label='RICH')
		self.ax.plot([-self.z_rich, self.z_rich], [self.r_out_rich, self.r_out_rich], color='lightsteelblue')
		self.ax.plot([self.z_rich, self.z_rich], [self.r_in_rich, self.r_out_rich], color='lightsteelblue')
		self.ax.plot([-self.z_rich, -self.z_rich], [self.r_in_rich, self.r_out_rich], color='lightsteelblue')

		if(self.color_areas):
			self.ax.fill_between([-self.z_rich, self.z_rich], [self.r_in_rich, self.r_in_rich], [self.r_out_rich, self.r_out_rich], color='lightsteelblue')

	def draw_FRICH(self):
		# A-side
		self.ax.plot([self.z_frich, self.z_frich], [self.r_in_frich, self.r_out_frich], color='tab:gray', label='fRICH')
		self.ax.plot([self.z_frich, self.z_frich + self.length_frich], [self.r_in_frich, self.r_in_frich], color='tab:gray')
		self.ax.plot([self.z_frich, self.z_frich + self.length_frich], [self.r_out_frich, self.r_out_frich], color='tab:gray')
		self.ax.plot([self.z_frich + self.length_frich, self.z_frich + self.length_frich], [self.r_in_frich, self.r_out_frich], color='tab:gray')

		# C-side
		self.ax.plot([-self.z_frich, -self.z_frich], [self.r_in_frich, self.r_out_frich], color='tab:gray')
		self.ax.plot([-self.z_frich, -self.z_frich - self.length_frich], [self.r_in_frich, self.r_in_frich], color='tab:gray')
		self.ax.plot([-self.z_frich, -self.z_frich - self.length_frich], [self.r_out_frich, self.r_out_frich], color='tab:gray')
		self.ax.plot([-self.z_frich - self.length_frich,-self.z_frich - self.length_frich], [self.r_in_frich, self.r_out_frich], color='tab:gray')

		if(self.color_areas):
			self.ax.fill_between([self.z_frich, self.z_frich + self.length_frich], [self.r_in_frich, self.r_in_frich], [self.r_out_frich, self.r_out_frich], color='tab:gray')
			self.ax.fill_between([-self.z_frich, -self.z_frich - self.length_frich], [self.r_in_frich, self.r_in_frich], [self.r_out_frich, self.r_out_frich], color='tab:gray')

	def draw_FCT_RICH(self):
		# Draw gas container
		self.ax.plot([self.z_fct_rich, self.z_fct_rich + self.fct_rich_radiator_length], [self.r_in_fct_rich, self.r_in_fct_rich], color='lime', label='FCT RICH')
		self.ax.plot([self.z_fct_rich, self.z_fct_rich + self.fct_rich_radiator_length], [self.r_out_fct_rich_lowz, self.r_out_fct_rich_highz], color='lime')
		self.ax.plot([self.z_fct_rich, self.z_fct_rich], [self.r_in_fct_rich, self.r_out_fct_rich_lowz], color='lime')
		self.ax.plot([self.z_fct_rich + self.fct_rich_radiator_length, self.z_fct_rich + self.fct_rich_radiator_length], [self.r_in_fct_rich, self.r_out_fct_rich_highz], color='lime')
		# Draw service container

		if(self.color_areas):
			self.ax.fill_between([self.z_fct_rich, self.z_fct_rich + self.fct_rich_radiator_length], [self.r_in_fct_rich, self.r_in_fct_rich], [self.r_out_fct_rich_lowz, self.r_out_fct_rich_highz], color='lime')
	
	def draw_FOCAL(self):
		self.ax.plot([self.z_focal, self.z_focal], [self.r_in_lowz_focal, self.r_out_lowz_focal], color='gold', label='FOCAL')
		self.ax.plot([self.z_focal, self.z_focal + self.length_focal], [self.r_in_lowz_focal, self.r_in_highz_focal], color='gold')
		self.ax.plot([self.z_focal, self.z_focal + self.length_focal], [self.r_out_lowz_focal, self.r_out_highz_focal], color='gold')
		self.ax.plot([self.z_focal + self.length_focal, self.z_focal + self.length_focal], [self.r_in_highz_focal, self.r_out_highz_focal], color='gold')

		if(self.color_areas):
			self.ax.fill_between([self.z_focal, self.z_focal + self.length_focal], [self.r_in_lowz_focal, self.r_in_highz_focal], [self.r_out_lowz_focal, self.r_out_highz_focal], color='gold')

	def draw_passive_material(self):
		# Draw iris vacuum vessel - inner tube
		self.ax.plot([-self.z_its[2] - 0.02, self.z_its[2] + 0.02], [self.r_trk[0] - 0.0002, self.r_trk[0] - 0.0002], color='black', label='Iris VacV')
		# Draw iris vacuum vessel - outer tube
		self.ax.plot([-self.z_its[2] - 0.02, self.z_its[2] + 0.02], [self.r_trk[2] + 0.002, self.r_trk[2] + 0.002], color='black')
		# Draw iris vacuum vessel - neg z side wall
		self.ax.plot([-self.z_its[2] - 0.02, -self.z_its[2] - 0.02], [self.r_trk[0] - 0.0002, self.r_trk[2] + 0.002], color='black')
		# Draw iris vacuum vessel - pos z side wall
		self.ax.plot([self.z_its[2] + 0.02, self.z_its[2] + 0.02], [self.r_trk[0] - 0.0002, self.r_trk[2] + 0.002], color='black')

		# Draw horizontal services for the vacuum vessel
		self.ax.plot([self.z_midpoint_services_vacv - self.z_services_vacv, self.z_midpoint_services_vacv + self.z_services_vacv], [self.r_in_services_vacv, self.r_in_services_vacv], color = 'mediumslateblue', label='Iris services')
		self.ax.plot([self.z_midpoint_services_vacv - self.z_services_vacv, self.z_midpoint_services_vacv + self.z_services_vacv], [self.r_out_services_vacv, self.r_out_services_vacv], color = 'mediumslateblue')
		self.ax.plot([self.z_midpoint_services_vacv - self.z_services_vacv, self.z_midpoint_services_vacv - self.z_services_vacv], [self.r_in_services_vacv, self.r_out_services_vacv], color = 'mediumslateblue')
		self.ax.plot([self.z_midpoint_services_vacv + self.z_services_vacv, self.z_midpoint_services_vacv + self.z_services_vacv], [self.r_in_services_vacv, self.r_out_services_vacv], color = 'mediumslateblue')

		if(self.color_areas):
			self.ax.fill_between([self.z_midpoint_services_vacv - self.z_services_vacv, self.z_midpoint_services_vacv + self.z_services_vacv], [self.r_in_services_vacv, self.r_in_services_vacv], [self.r_out_services_vacv, self.r_out_services_vacv], color='mediumslateblue')

		# Draw vertical services for the vacuum vessel
		self.ax.plot([self.z_its[2] + 0.03, self.z_its[2] + 0.03], [self.r_trk[0] - 0.0002, self.r_in_services_vacv], color = 'mediumslateblue')
		self.ax.plot([self.z_vacuumVesselWall - 0.02, self.z_vacuumVesselWall - 0.02], [self.r_trk[0] - 0.0002, self.r_in_services_vacv], color = 'mediumslateblue')
		self.ax.plot([self.z_its[2] + 0.03, self.z_vacuumVesselWall - 0.02], [self.r_trk[0] - 0.0002, self.r_trk[0] - 0.0002], color = 'mediumslateblue')
		self.ax.plot([self.z_its[2] + 0.03, self.z_vacuumVesselWall - 0.02], [self.r_in_services_vacv, self.r_in_services_vacv], color = 'mediumslateblue')

		if(self.color_areas):
			self.ax.fill_between([self.z_its[2] + 0.03, self.z_vacuumVesselWall - 0.02], [self.r_trk[0] - 0.0002, self.r_trk[0] - 0.0002], [self.r_in_services_vacv, self.r_in_services_vacv], color='mediumslateblue')

	def draw_coldplate(self):
		self.ax.plot([-self.z_coldplate, self.z_coldplate], [self.r_coldplate, self.r_coldplate], color='limegreen', label='Cold plate')

	def draw_pseu_lines(self):
		z_origin_list = []
		y_origin_list = []
		z_disp_list = []
		y_disp_list = []
		labelp = r'$\eta$ = '
		for z_pos in self.pseurap_origins:
			for eta in self.pseudorapidity_list:
				ang = self.pseurap_to_angle(abs(eta))
				z = []
				y = []
				if(eta < 0):
					z = [-self.xright, z_pos]
					y = [(self.xright + z_pos) * math.tan(ang), 0]
				if(eta > 0):
					z = [z_pos, self.xright]
					y = [0, (self.xright - z_pos) * math.tan(ang)]
				if(z_pos == 0.):
					z_origin_list.append(z)
					y_origin_list.append(y)
					labelp += f'{str(eta)}, '
				else:
					z_disp_list.append(z)
					y_disp_list.append(y)

					
		for i in range(0, len(z_origin_list)):
			if(i == 0):
				self.ax.plot(z_origin_list[i], y_origin_list[i], 'r--', label=labelp[:-2])
			else:
				self.ax.plot(z_origin_list[i], y_origin_list[i], 'r--')

		for i in range(0, len(z_disp_list)):
			self.ax.plot(z_disp_list[i], y_disp_list[i], 'k-')

toymodel = Toy_FCT()
# toymodel.set_inner_disc_radii(4, 1)
toymodel.print_disc_details(0)
toymodel.print_disc_details(1)
# toymodel.fct_rich_in_front_of_fct(0.7)
# toymodel.write_layout(0, "FCT") # Writes layout that can be used for O2
# toymodel.window = True
# toymodel.open_iris()
toymodel.draw_all()

toymodel.grid()
toymodel.xleft = 0
# toymodel.xright = 3
# toymodel.yup = 0.4
toymodel.show()
