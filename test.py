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