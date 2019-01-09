_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
#import BEG
import clusters
import js
import m_structure
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from sklearn.linear_model import Ridge
from sklearn import linear_model

def read_Js(num_Js):
    js_file = open('output','r')
    Js = []
    lines = js_file.readlines()
    for i in range(num_Js):
        line = lines[i+1]
        line = line.split()
        Js.append(float(line[2]))
    return Js

def ridge_regression(data, predictors, alpha, models_to_plot={}):
    # Fit the model
    ridgereg = linear_model.Ridge(alpha=alpha,normalize=True,fit_intercept=True)
    ridgereg.fit(predictors,data)
    y_pred = ridgereg.predict(predictors)

    # Check if a plot is to be made for the entered alpha
    if alpha in models_to_plot:
        plt.subplot(models_to_plot[alpha])
        plt.tight_layout()
        plt.plot(data,y_pred,'.')
        #plt.plot(data['x'],data['y'],'.')
        plt.xlim(min(min(data),min(y_pred))-5,max(max(data),max(y_pred))+5)
        plt.ylim(min(min(data),min(y_pred))-5,max(max(data),max(y_pred))+5)
        plt.title('alpha: %.3g'%alpha)

    # Return result in predefined format
    rss = sum(np.power(y_pred-data, 2))
    #print(rss[0,0])
    ret = [ rss[0,0] ]
    ret.extend([ridgereg.intercept_])
    ret.extend(ridgereg.coef_[0])
    return ret

# Calculates the fitting parameters using a ridge regression with autocorrelation
def ridge_simple(m_structure_list,alpha,Cluster_rules,J_rules):
    predictors = []         # used to be a
    data = []               # used to be b
    for i in range(len(m_structure_list)):
        mat = m_structure_list[i]
        
        ## find indices for mag 3NN J-terms and set those sums to zero here.
        #for j in range(len(J_rules)):
            #if ("mag" in J_rules[j].name) and ("3NN" in J_rules[j].name):
                #mat.J_sums[j]=0
        
        row = mat.Cluster_sums + mat.J_sums
        for j in range(len(row)):
            row[j] *= mat.weight
        predictors.append(row)
        data.append(mat.enrg * mat.weight)
    predictors = np.matrix(predictors)
    data = np.transpose(np.matrix(data))

# evaluate coefficients for different regularization parameters and plot results
# would like to make this label the different curves with the parameter name
    n_alphas = 510
    alpha_list = np.logspace(-15,10,n_alphas)
    
    # initialize dataframe for storing coefficients
    col = ['rss','intercept']+['%s'%Cluster_rules[i].name.rstrip() for i in range(len(Cluster_rules))]+['%s'%J_rules[i].name.rstrip() for i in range(len(J_rules))]
    ind = ['alpha_%.2g'%alpha_list[i] for i in range(len(alpha_list))]
    coef_matrix_ridge = pd.DataFrame(index=ind,columns=col)
    
    plt.figure
    models_to_plot = {alpha_list[0]:231,alpha_list[10]:232,alpha_list[20]:233,alpha_list[30]:234,alpha_list[40]:235,alpha_list[50]:236}
    for i in range(len(alpha_list)):
        coef_matrix_ridge.iloc[i,] = ridge_regression(data,predictors,alpha_list[i],models_to_plot)
    plt.savefig('testing_regularization.pdf')

    pd.options.display.float_format = '{:,.2g}'.format
    f1=open('RegularizationParameters.txt', 'w')
    print(coef_matrix_ridge, file=f1)
    f1.close()

#    coefs = []
#    scores = []
#    for j in alpha_list:
#        ridge_fit = linear_model.Ridge(alpha=j,fit_intercept=False,normalize=True)
#        ridge_fit.fit(a,b)
#        coefs.append(ridge_fit.coef_[0])
#        score=ridge_fit.score(a,b)
#        scores.append(score)
    plt.figure(2)
    ax = plt.gca()
    c_plot=[]
    for k in coef_matrix_ridge.keys():
        if (('NN' in k) or ('mono' in k)) and np.max(np.abs(coef_matrix_ridge[k])) > 1e-8:
           c_plot.append(k)
    #c_plot.append('intercept')
    for k in c_plot:
        ax.plot(alpha_list, coef_matrix_ridge[k],label=k)
    ax.set_xscale('log')
    ax.set_xlim(ax.get_xlim()[::-1])  # reverse axis
    ax.legend(loc=(1.05,0.))
    plt.xlabel('regularization parameter alpha')
    plt.ylabel('fitted coefficients')
    plt.title('Ridge coefficients as a function of regularization')
    plt.axis('tight')
    plt.savefig('regularization.pdf',bbox_inches='tight')

    plt.figure(3)
    bx = plt.gca()
    #bx.plot(alpha_list, scores)
    bx.plot(alpha_list,coef_matrix_ridge['rss'])
    bx.set_xscale('log')
    bx.set_yscale('log')
    bx.set_xlim(bx.get_xlim()[::-1])  # reverse axis
    plt.xlabel('regularization parameter alpha')
    plt.ylabel('Residual Sum of Squares')
    plt.title('RSS as a function of regularization')
    plt.axis('tight')
    plt.savefig('scores.pdf')
    
    # determine regularization parameter using cross-validation
    ridge_fit = linear_model.RidgeCV(alphas=alpha_list,fit_intercept=False,normalize=True)
    ridge_fit.fit(predictors,data)
    print('Selected regularization parameter using cross validation: ',ridge_fit.alpha_)

    #ridge_fit = linear_model.Ridge(alpha=0.00001,fit_intercept=False,normalize=True) ###################################
    #ridge_fit.fit(predictors,data)

    #ridge_fit = linear_model.Ridge(alpha=1e-10,fit_intercept=False)
    #ridge_fit.fit(a,b)

    Js = ridge_fit.coef_[0]
    #intercept = ridge_fit.intercept_[0]
    intercept = 0
    #predict=ridge_fit.predict(a)
    #print(predict)
    #JS_list = []
    #Js_0 = Js[0]
    #for i in range(len(Js_0)):
    #Js_0
    #JS_list.append(Js_0[i])
    #print(Js)
    #print(ridge_fit.predict(a))
    #print(ridge_fit.intercept_)
    return Js,intercept







def write_fitting_parameters(structures, clusters_list, j_list, Js, intercept, limit):
    label = []
    #    for i in range(len(beg_list)):
    #label.append(beg_list[i].name)
    for i in range(len(clusters_list)):
        label.append(clusters_list[i].name)
    for i in range(len(j_list)):
        label.append(j_list[i].name)
    path = 'FittingParameters'
    file = open(path, 'w')
    file.write("Fitting Parameters\n")
    for i in range(len(Js)):
        line = label[i].strip() + " = " + str(Js[i]) + "\n"
        line = line.replace("[", "")
        line = line.replace("]", "")
        file.write(line)
    file.write("\n")
    file.write("Name".ljust(20)+"Actual Energy".ljust(20)+"Predicted Energy".ljust(20)+"Magnetic Order".ljust(20)+"Phase".ljust(20)+"Weight".ljust(20)+"\n")
    for i in range(len(structures)):         # seems like energy evaluation should just be a function we can call
        mat = structures[i]
        #        if mat.phase_name != "pmmmm":
        file.write(str(mat.name).ljust(20) + str(mat.enrg).ljust(20))
        new_enrg = intercept
        #for j in range(len(beg_list)):
        #new_enrg += Js[j] * structures[i].BEG_sums[j]
        for k in range(len(clusters_list)):
            new_enrg += Js[k] * structures[i].Cluster_sums[k]
        for l in range(len(j_list)):
            new_enrg += Js[len(clusters_list) + l] * structures[i].J_sums[l]
        line = str(new_enrg)
        line = line.replace("[", "")
        line = line.replace("]", "")
        file.write(line.ljust(20) + mat.mag_phase.ljust(20) + mat.phase_name.ljust(20) + str(mat.weight).ljust(20) + "\n")
    file.close()






# plotting function not written very generally here, just does NiMnIn
def plot_data3(M_structures, clusters_list, j_list, Js, intercept, limit,hull):
    plt.figure(4)           # do we really need to keep count of the number of figures?
    colors = {'mart': 'g', 'aus': 'b', 'pre-mart': 'r'}
    markers = {'none': 'x', 'FM': 'o', 'AFM': 's', 'spin disordered':'D'}
    labels = ['NiMn', 'Ni4Mn3In1', 'Ni2MnIn']
    actual_labels = ['Ni$_8$Mn$_8$', 'Ni$_8$Mn$_6$In$_2$', 'Ni$_8$Mn$_4$In$_4$']
    e_comp0 = []
    e_comp50 = []
    for i in range(len(M_structures)):
        if M_structures[i].phase_name != "pmmmm":
            structure = M_structures[i]
            comp = structure.composition[2]/structure.composition[0]
            if comp == 0:
                e_comp0.append(structure.enrg)
            if comp == .5:
                e_comp50.append(structure.enrg)
    comp0_min = min(e_comp0)
    comp50_min = min(e_comp50)
    offset = (comp50_min-comp0_min)/.5
    #offset = 0
    enrg_list = []
    for i in range(len(M_structures)):
        if M_structures[i].phase_name != "pmmmm":
            comp = M_structures[i].composition[2]/M_structures[i].composition[0]
            if hull == True:
                enrg_list.append(M_structures[i].enrg)
            else:
                enrg_list.append(M_structures[i].enrg-(offset*comp + comp0_min))
    new_enrg_list = []
    for i in range(len(M_structures)):
        mat = M_structures[i]
        comp = mat.composition[2]/mat.composition[0]
        if mat.phase_name != "pmmmm":
            new_enrg = intercept
            #for j in range(len(beg_list)):
            #new_enrg += Js[j] * M_structures[i].BEG_sums[j]        # ugh - here evaluating it again, now above hull
            for k in range(len(clusters_list)):
                new_enrg += Js[k] * M_structures[i].Cluster_sums[k]
            for l in range(len(j_list)):
                new_enrg += Js[len(clusters_list) + l] * M_structures[i].J_sums[l]
            new_enrg -= offset*comp + comp0_min                         #
            new_enrg_list.append(new_enrg)
    x = 0
    x2_itter = -.5
    x4_itter = -.5
    x6_itter = -.5
    for i in range(len(enrg_list)):
        if M_structures[i].phase_name != "pmmmm":
            if int(M_structures[i].composition[1]) == 8:
                x = 1.5
                x2_itter += .08
                x_itter = x2_itter
            if int(M_structures[i].composition[1]) == 6:
                x = 4.5
                x4_itter += .08
                x_itter = x4_itter
            if int(M_structures[i].composition[1]) == 4:
                x = 9.5
                x6_itter += .08
                x_itter = x6_itter
            if M_structures[i].phase_name == "aust":
                c = 'b'
            if M_structures[i].phase_name == "mart":
                c = 'g'
            if M_structures[i].phase_name == "pm":
                c = 'r'
            if M_structures[i].mag_phase == "fm":
                m = 'o'
            if M_structures[i].mag_phase == "afm":
                m = 's'
            if M_structures[i].mag_phase == "pera":
                m = '^'
            if M_structures[i].mag_phase == "sd":
                m = 'D'
            if M_structures[i].mag_phase == "NA":
                m = 'x'
            y = float(enrg_list[i])
            if M_structures[i].phase_name != 'pm':
                #plt.plot([x + x_itter, x + x_itter], [y, float(new_enrg_list[i])], lw=1, color="k")
                plt.plot(x + x_itter, y, lw=0, markersize=8, marker=m, color=c)
                plt.plot(x + x_itter, float(new_enrg_list[i]), lw=0, markersize=8, marker=".", color="r")
    for c in colors.keys():
        for m in markers.keys():
            plt.plot([],[],label=c+', '+m,lw=0,markersize=12,marker=markers[m],color=colors[c])
    #plt.xlim(0,8)
    #plt.ylim(-1,9)
    plt.rc('lines', linewidth=1)
    plt.title("NiMn -- Ni$_2$MnIn Composition Energies", fontsize=16)
    plt.xlabel("Composition", fontsize=14)
    plt.ylabel("Energy above Hull (eV/16 atoms)", fontsize=14)
    #plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)
    plt.xticks([2,5.5,10],actual_labels, rotation='horizontal',fontsize=14)
    plt.savefig('FittedDataSummary.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
