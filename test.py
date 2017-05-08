#
# fig = plt.figure(5)
#                 ax = fig.add_subplot(111, projection='3d')
#                 xs = []
#                 ys = []
#                 zs = []
#                 cs = []
#                 ms = []
#                 us = []
#                 vs = []
#                 ws = []
#                 ws2 = []
#                 for i in range(supercell_obj.i_length):
#                     for j in range(supercell_obj.j_length):
#                         for k in range(supercell_obj.k_length):
#                             if np.mod(k,2) == 0:
#                                 offset = 0
#                             else:
#                                 offset = .5
#                             site = [i,j,k]
#                             pos = supercell_obj.get_site_pos(site)
#                             xs.append(pos[0]+offset)
#                             ys.append(pos[1]+offset)
#                             zs.append(pos[2]*.5)
#                             us.append(0)
#                             vs.append(0)
#                             ws2.append(supercell_obj.get_site_spin(site))
#                             ws.append(supercell_obj.get_site_phase(site))
#                             if supercell_obj.get_site_species(site) == 0:
#                                 cs.append('g')
#                             if supercell_obj.get_site_species(site) == 1:
#                                 cs.append('r')
#                             if supercell_obj.get_site_species(site) == 2:
#                                 cs.append('b')
#                             if supercell_obj.get_site_pos(site) in cluster:
#                                 ms.append('s')
#                             else: ms.append('o')
#                 for i in range(supercell_obj.num_sites):
#                     ax.quiver(xs[i],ys[i],zs[i],us[i],vs[i],ws[i],pivot='middle',length=.5,color='b')
#                     ax.quiver(xs[i],ys[i],zs[i],us[i],vs[i],ws2[i],pivot='middle',length=.25,color='c')
#                     ax.scatter(xs[i],ys[i],zs[i],c=cs[i],marker=ms[i],s = 50)
#                 plt.show()






# def eval_site_new(site,supercell_obj,Cluster_rules,J_ruels,Js,T,B_ext):
#     Kb = .000086173324
#     uB = .000057883818012 #5.7883818012(26)×10−5 eV/Tesla
#     g = 2
#     total_Ham = 0
#     site_phase = supercell_obj.get_site_phase(site)
#     J,K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
#     for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
#         if supercell_obj.get_neighbor_order(site,neighbor) == 1:
#             neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
#             total_Ham += J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2)
#     total_Ham += Kb*T*np.log(8)*(site_phase**2)
#     total_Ham += -g*uB*B_ext*supercell_obj.get_site_spin(site)*1
#     return total_Ham

#
# def eval_cluster(supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_ruels,Js,T):
#     Kb = .000086173324
#     total_H = 0
#     total_H2 = 0
#     total_H3 = 0
#     total_H4 = 0
#     count = 0
#     count2 = 0
#     count3 = 0
#     for i in range(len(links)):
#         site = links[i]
#         site_phase = supercell_obj.get_site_phase(site)
#         J,K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
#         for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
#             if supercell_obj.get_neighbor_order(site,neighbor) == 1:
#                 if supercell_obj.get_neighbor_pos(site,neighbor) in links:
#                     count += 1
#                     neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
#                     total_H += J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2)
#         total_H += Kb*T*np.log(8)*(site_phase**2)
#
#         for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
#             if supercell_obj.get_neighbor_order(site,neighbor) == 1:
#                 count2 += 1
#                 neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
#                 total_H2 += J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2)
#         total_H2 += Kb*T*np.log(8)*(site_phase**2)
#
#         for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
#             if supercell_obj.get_neighbor_order(site,neighbor) == 1:
#                 if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
#                     count3 += 1
#                     neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
#                     total_H3 += J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2)
#
#         total_H4 += Kb*T*np.log(8)*(site_phase**2)
#
#     print([[total_H,total_H2,total_H3],[count,count2,count3],[total_H4]])
#     return total_H