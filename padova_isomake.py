import numpy as np
import scipy as sp
import scipy.ndimage as snd
import scipy.interpolate as si
import gzip as gz
import glob as gb




def iso_interp(filename, metallicity, metal_weight, output_obj, bands_dict):
    print metallicity, metal_weight

    logT_min=3.3
    logT_max=4.6
    logT_step=0.01
    
    logg_min=-1.5
    logg_max=5.5
    logg_step=0.01

    logT_edges=np.arange(logT_min-logT_step/2., logT_max+logT_step,logT_step)
    logg_edges= np.arange(logg_min-logg_step/2., logg_max+logg_step,logg_step)

    interp_points_grid=np.mgrid[logT_min:logT_max+logT_step:logT_step, logg_min:logg_max+logg_step:logg_step]
    interp_points=np.array([interp_points_grid[0].flatten(), interp_points_grid[1].flatten()])
    
    # Find columns required
    
    with open(filename, 'r') as f:
        for x in (100)
            header_line=f.readline().split()
            if "M_ini" in header_line:
                break    
    
    Mi_col = header_line.index("M_ini")-1
    logage_col = header_line.index("log(age/yr)")-1
    logTe_col = header_line.index("logTe")-1
    logg_col = header_line.index("logG")-1
    
    photom_cols={}
    
    for band in bands_dict.keys():
        if band in header_line:
            photom_cols[bands_dict[band] ] =  = header_line.index(band)-1
        

    # Read in data - need Z, age, Mi, logT, logg, r, i, Ha

    iso_data=np.loadtxt( filename, usecols=[Mi_col, logage_col, logTe_col, logg_col].extend(photom_cols.values() ) )

    # Set metallicity

    iso_data[:,0]=metallicity

    # work out point weight

    weights=np.zeros([iso_data.shape[0],1])
    for i in range(iso_data.shape[0]):
        if i==0:
            weights[i]=(iso_data[i+1,2]-iso_data[i,2])*np.power(10,iso_data[i,1])
        elif i==iso_data.shape[0]-1:
            weights[i]=(iso_data[i,2]-iso_data[i-1,2])*np.power(10,iso_data[i,1])
        elif iso_data[i,2]>iso_data[i+1,2]:
            weights[i]=(iso_data[i,2]-iso_data[i-1,2])*np.power(10,iso_data[i,1])
        elif iso_data[i-1,2]>iso_data[i,2]:
            weights[i]=(iso_data[i+1,2]-iso_data[i,2])*np.power(10,iso_data[i,1])
        else:
            weights[i]=(iso_data[i+1,2]-iso_data[i-1,2])*np.power(10,iso_data[i,1])
    weights1=weights[:].flatten()
    weights.reshape(iso_data.shape[0],1)
    iso_data=np.append(iso_data, weights, axis=1)

    # Index data

    binned_data=[ [ [] for i in range(13) ] for i in range(13) ]

    for i in range(iso_data.shape[0]):

        try:
            binned_data[int(np.floor((iso_data[i,3]-3.4)/0.1))][int(np.floor((iso_data[i,4]+0.5)/0.5))].append(iso_data[:])
        except IndexError:
            print datum[3], datum[4] 



    # Re grid in hi-res

    metals=np.zeros(interp_points[0].shape)+metallicity

#    logage_spline=si.SmoothBivariateSpline(iso_data[:,3], iso_data[:,4], iso_data[:,1], w=weights)
#    #logage_rbf=si.Rbf(iso_data[:,3], iso_data[:,4], iso_data[:,1])
#    for i in range(0, interp_points[0].size, 10):
#        logage=logage_spline(interp_points[0,i:i+10], interp_points[1,i:i+10])
##        logage=logage_rbf(interp_points[0,i:i+250], interp_points[1,iIi+250])

#    logage=si.griddata(iso_data[:,3:5], iso_data[:,1], interp_points.T, method='cubic', fill_value=-99)
#    Mi=si.griddata(iso_data[:,3:5], iso_data[:,2], interp_points.T, method='cubic', fill_value=-99)
#    
#    r=si.griddata(iso_data[:,3:5], iso_data[:,5], interp_points.T, method='cubic', fill_value=-99)
#    i=si.griddata(iso_data[:,3:5], iso_data[:,6], interp_points.T, method='cubic', fill_value=-99)
#    ha=si.griddata(iso_data[:,3:5], iso_data[:,7], interp_points.T, method='cubic', fill_value=-99)
    
    logage=np.zeros(interp_points.T.shape[0])
    Mi=np.zeros(interp_points.T.shape[0])
    r=np.zeros(interp_points.T.shape[0])
    i=np.zeros(interp_points.T.shape[0])
    ha=np.zeros(interp_points.T.shape[0])

    for it, point in enumerate(interp_points.T):
#        print point, int(np.floor((point[0]-3.4)/0.1)), int(np.floor((point[1]+0.5)/0.5))
#        iso_list=binned_data[int(np.floor((point[0]-3.4)/0.1))][int(np.floor((point[1]+0.5)/0.5))]

        shortlist=iso_data[np.power(iso_data[:,3]-point[0],2)*4. + np.power(iso_data[:,4]-point[1],2)< 0.01]
        if shortlist.size==0:
            logage[it]=-99
            Mi[it]=-99
            r[it]=-99
            i[it]=-99
            ha[it]=-99
        else:
            short_weights1=shortlist[:,-1]/(1E-9+np.power(shortlist[:,3]-point[0],2.)*36.+np.power(shortlist[:,4]-point[1],2.))
            short_weights_sum1=np.sum(short_weights1)
            short_weights2=1/(1E-9+np.power(shortlist[:,3]-point[0],2.)*36.+np.power(shortlist[:,4]-point[1],2.))
            short_weights_sum2=np.sum(short_weights2)

            logage[it]=np.sum(shortlist[:,1]*short_weights1)/short_weights_sum1
            Mi[it]=np.sum(shortlist[:,2]*short_weights1)/short_weights_sum1
            r[it]=np.sum(shortlist[:,5]*short_weights2)/short_weights_sum2
            i[it]=np.sum(shortlist[:,6]*short_weights2)/short_weights_sum2
            ha[it]=np.sum(shortlist[:,7]*short_weights2)/short_weights_sum2
    
    inner_counts, xe, ye=np.histogram2d(iso_data[:,3], iso_data[:,4], bins=[logT_edges, logg_edges])
    outer_counts=(snd.filters.uniform_filter(inner_counts, size=(2,5))*25+0.1).astype(int)
    inner_Jac, xe, ye=np.histogram2d(iso_data[:,3], iso_data[:,4], bins=[logT_edges, logg_edges], weights=weights1)
    outer_Jac=(snd.filters.uniform_filter(inner_Jac, size=(2,5))*25)

    outer_Jac[outer_Jac<0.01]=0.
    
    # Save to file

    np.savetxt(output_obj, np.array([metals, Mi, logage, interp_points[0], interp_points[1], outer_Jac.flatten(), r,i,ha, inner_counts.flatten(), outer_counts.flatten() ]).T, fmt=['%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3e', '%.3f', '%.3f', '%.3f', '%i', '%i'])
    
    
    
    
    
    
    
def padova_interpolated_isomake(folders, bands, output_file):
