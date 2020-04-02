

rleak1 = []
rleak2 = []

gkleak_array1=[]
gkleak_array2=[]

gs_array1=['1.0']

t_array1=['10']

gh_array1=[]
gh_array2=[]

#Use these values og gh and and potassium leak conductance for Figure8
#P = [(28,80),(29,82.5),(30,85),(31,87.5),(32,90),(33,92.5),(34,95),(35,97.5),(36,101),(37,103),(38,105),(39,107),(40,109),(41,111.25),(43,113.5),(43,115.75),(44,118)]

for g in range(40,180,1):
	gkleak_array1.append(format(g/10000.000000000,'.6f'))
	gkleak_array2.append(str(g))
	# k_RE = 800.00*(1 - (g - 100.00)/100.00)

	# rleak1.append(format(k_RE/10000.000000000,'.6f'))
	# rleak2.append(str(g))


for g in range(8,50,1):
	gh_array1.append(format(g/100.000000000,'.6f'))
	gh_array2.append(str(g))


calciumshift1 = []
for c in range(100,101,1): #chnage range from 50 to 150 (%) for Figure8
	calciumshift1.append(format(c/100.000000000,'.2f'))

calciumshift2 = []
for c in range(100,101,1):
	calciumshift2.append(format(c/100.000000000,'.2f'))


for c1 in calciumshift1:
	for c2 in calciumshift2:
		for c in ['0.00024']:
			for l in range(0,len(gkleak_array2)):

				for t in range(0,len(t_array1)):

					for g in range(0,len(gh_array1)):

						NSETS_CONF="./nsets_network_mAch_Ca"+c+"_gh"+gh_array2[g]+"_leak"+gkleak_array2[l]+"_tepsp10_calciumshift1"+c1+"_calciumshift2"+c2+"_RE.isf"
						strtofile = ""


						for i in range(0,2):

							NEURON="\" HTC_cell_"+str(i)+" \" \n dxdt:9 ,"
							NEURON+="m_Na_TC:0.2,"
							NEURON+="h_Na_TC:0.6,"
							NEURON+="n_K_TC:0.4,"
							NEURON+="h_t_TC:0.024,"
							NEURON+="m_ahp_HTC:0.05,"
							NEURON+="r_H_HTC:0.5,"
							NEURON+="h_tht_HTC:0.40,"
							NEURON+="v:-60.0,"
							NEURON+="Ca_conc:0.00024,"
							NEURON+="m_Na_TC_stoch:-1,"
							NEURON+="I_Na_RE:0,"
							NEURON+="I_K_RE:0,"
							NEURON+="Ca_slow_switch:0,"
							NEURON+="I_T_RE:0,"
							NEURON+="v_mag:0.01,"
							NEURON+="I_AHP_HTC:0,"
							NEURON+="m_Na_RE:0.0,"
							NEURON+="n_K_RE:0.0,"
							NEURON+="h_Na_RE:0.0,"
							NEURON+="m_t_RE:0.0,"
							NEURON+="h_t_RE:0.0,"
							NEURON+="I_T_TC:0,"
							NEURON+="eq_Ca_fast:"+c+","
							NEURON+="eq_Ca_slow:0.00024,"
							NEURON+="tau_fast:3,"
							NEURON+="tau_slow:9999999999999999999999,"
							NEURON+="gkleak:"+gkleak_array1[l]+","
							NEURON+="h_Na_TC_stoch:-1,"
							NEURON+="n_K_TC_stoch:-1,"
							NEURON+="eq_Ca:0.00024,"
							NEURON+="tau_Ca:99999999,"
							NEURON+="h_t_TC_stoch:-1,"
							NEURON+="I_THT_HTC:0,"
							NEURON+="Ca_conc_stoch:-1,"
							NEURON+="calcium_shift1:"+c1+","
							NEURON+="calcium_shift2:"+c2+","
							NEURON+="m_ahp_HTC_stoch:-1,"
							NEURON+="r_H_HTC_stoch:-1,"
							NEURON+="h_tht_HTC_stoch:-1,"
							NEURON+="v_stoch:1,"
							NEURON+="I_Syn:0,"
							NEURON+="I_Syn_I:0,"
							NEURON+="I_Syn_E:0,"
							NEURON+="T:0,"
							NEURON+="last_spike:-999999999,"
							NEURON+="g_THT:12,"
							NEURON+="I_Syn_GJ:0,"
							NEURON+="gh:"+gh_array1[g]+","
							NEURON+="g_gj:0.005,"
							NEURON+="iext:0; \n"
							strtofile+=NEURON



						for i in range(2,10):
							NEURON="\" TC_cell_"+str(i)+" \" \n	dxdt:9,"
							NEURON+="m_Na_TC:0.1,"
							NEURON+="n_K_TC:0.4,"
							NEURON+="h_Na_TC:0.6,"
							NEURON+="h_t_TC:0.2,"
							NEURON+="o_h_TC:0.2,"
							NEURON+="c_h_TC:0.8,"
							NEURON+="Ca_conc:0.00024,"
							NEURON+="p_h_TC:0.5,"
							NEURON+="v:-56,"
							NEURON+="I_T_TC:0,"
							NEURON+="I_Syn:0,"
							NEURON+="I_Syn_I:0,"
							NEURON+="I_Syn_E:0,"
							NEURON+="I_Na_RE:0,"
							NEURON+="I_K_RE:0,"
							NEURON+="I_T_RE:0,"
							NEURON+="T_next_epsp:4,"
							NEURON+="m_Na_RE:0.0,"
							NEURON+="n_K_RE:0.0,"
							NEURON+="h_Na_RE:0.0,"
							NEURON+="m_t_RE:0.0,"
							NEURON+="h_t_RE:0.0,"
							NEURON+="v_mag:1,"
							NEURON+="T_last_epsp:-999999999.0,"
							NEURON+="T_last_ipsp:-999999999.0,"
							NEURON+="T_next_ipsp:99999999.0,"
							NEURON+="T_next_epsp:100,"
							NEURON+="Ca_slow_switch:0,"
							NEURON+="eq_Ca_fast:"+c+","
							NEURON+="t_ipsp:9999999,"
							NEURON+="t_epsp:10,"
							NEURON+="eq_Ca_slow:0.00024,"
							NEURON+="tau_fast:20,"
							NEURON+="tau_slow:10,"
							NEURON+="calcium_shift1:"+c1+","
							NEURON+="calcium_shift2:"+c2+","
							NEURON+="m_Na_TC_stoch:-1,"
							NEURON+="n_K_TC_stoch:-1,"
							NEURON+="h_Na_TC_stoch:-1,"
							NEURON+="h_t_TC_stoch:-1,"
							NEURON+="eq_Ca:0.00024,"
							NEURON+="tau_Ca:3,"
							NEURON+="r_H_HTC:0.7777777,"
							NEURON+="o_h_TC_stoch:-1,"
							NEURON+="c_h_TC_stoch:-1,"
							NEURON+="I_THT_HTC:0,"
							NEURON+="Ca_conc_stoch:-1,"
							NEURON+="p_h_TC_stoch:-1,"
							NEURON+="last_spike:-999999999,"
							NEURON+="v_stoch:1,"
							NEURON+="T:0,"
							NEURON+="gs:1.0,"
							NEURON+="iext:0.0; \n"
							strtofile+=NEURON





						for i in range(10,20):
							NEURON="\" RE_cell_"+str(i)+" \" \n dxdt:7,"
							NEURON+="m_Na_RE:0.1,"
							NEURON+="n_K_RE:0.4,"
							NEURON+="h_Na_RE:0.6,"
							NEURON+="Ca_conc:0.00024,"
							NEURON+="m_t_RE:0.5,"
							NEURON+="h_t_RE:0.5,"
							NEURON+="v:-60.0,"
							NEURON+="I_T_TC:0,"
							NEURON+="I_Na_RE:0,"
							NEURON+="I_K_RE:0,"
							NEURON+="I_Syn:0,"
							NEURON+="I_Syn_I:0,"
							NEURON+="I_Syn_E:0,"
							NEURON+="gs:1.0,"
							NEURON+="I_T_RE:0,"
							NEURON+="v_mag:0.1,"
							NEURON+="T_last_epsp:-999999999.0,"
							NEURON+="T_next_epsp:1.0,"
							NEURON+="T_last_ipsp:-999999999.0,"
							NEURON+="T_next_ipsp:3.0,"
							NEURON+="I_AHP_HTC:0,"
							NEURON+="m_Na_RE_stoch:-1,"
							NEURON+="t_ipsp:10,"
							NEURON+="t_epsp:10,"
							NEURON+="eq_Ca_fast:"+c+","
							NEURON+="eq_Ca_slow:0.00024,"
							NEURON+="r_H_HTC:0.7777777,"
							NEURON+="tau_fast:3,"
							NEURON+="tau_slow:3,"
							NEURON+="n_K_RE_stoch:-1,"
							NEURON+="Ca_conc_stoch:-1,"
							NEURON+="h_Na_RE_stoch:-1,"
							NEURON+="m_t_RE_stoch:-1,"
							NEURON+="h_t_RE_stoch:-1,"
							NEURON+="calcium_shift1:"+c1+","
							NEURON+="calcium_shift2:"+c2+","
							NEURON+="T:0,"
							NEURON+="last_spike:-999999999,"
							NEURON+="v_stoch:1,"
							NEURON+="gkleak:0.08,"
							NEURON+="iext:0.0; \n"
							strtofile+=NEURON


						outfile = open(NSETS_CONF,"w")

						outfile.write(strtofile)
						outfile.close()



print(rleak1)





