import glob as gb
import gzip as gz
import numpy as np
import scipy as sp
import scipy.ndimage as snd
import scipy.interpolate as si
from sys import stdout



def replacement_gradient(x):

    diff=np.diff(x)
    out=np.zeros(x.shape)
    out[:-1]+=diff
    out[1:]+=diff
    out[1:-1]/=2.
    
    return out


def iso_interp(filenames, metallicity, metal_weight, output_obj, bands_dict, bands_ordered):
    """ iso_interp(filenames, metallicity, metal_weight, 
                   output_obj, bands_dict, bands_ordered)

        Interpolate isochrones defined in (initial Mass, log(age))
        space onto a grid of points in (log(T_eff), log(g)).

        Parameters
        ----------
        filenames : list(string)
            the files that contain the isochrones
        metallicity : float
            the metallicity of the isochrone
        metal_weight : float
            a probability weight assigned to the metallicity
        output_obj : 
            a file output object to save the interpolated 
            isochrones to.
        bands_dict : dict
            A dictionary containing the photometric bands
            to be used.
        bands_ordered : list
            The bands in their desired order (can be arbitrary)

        Notes
        -----
        Set up for use with Padova isochrones. May or may not 
        work with other isochrone libraries.
    """

    print "[M/H]={0:.3f} , weight = {1:.3f}".format(metallicity, metal_weight)
    if isinstance(filenames, basestring):
        filenames=[filenames]

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
    
    iso_data=[]
    photom_data={}    
    
    for filename in filenames:
    
        with open(filename, 'r') as f:
            for x in range(100):
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
                photom_cols[bands_dict[band] ] = header_line.index(band)-1
            

        # Read in data - need Z, age, Mi, logT, logg, r, i, Ha

        iso_data.append(np.loadtxt( filename, usecols=(Mi_col, logage_col, logTe_col, logg_col) ))
        print iso_data[-1].size, iso_data[0].size, filename

        assert iso_data[-1].size==iso_data[0].size
        

        for band in photom_cols.keys():
            photom_data[band]=np.loadtxt( filename, usecols=[photom_cols[band]] )
            
    iso_data=iso_data[-1]            

    # Set metallicity


    # work out point weight

    weights=np.zeros([iso_data.shape[0],1])
    for i in range(iso_data.shape[0]):
        if i==0:
            weights[i]=(iso_data[i+1,0]-iso_data[i,0])*np.power(10,iso_data[i,1])
        elif i==iso_data.shape[0]-1:
            weights[i]=(iso_data[i,0]-iso_data[i-1,0])*np.power(10,iso_data[i,1])
        elif iso_data[i,0]>iso_data[i+1,0]:
            weights[i]=(iso_data[i,0]-iso_data[i-1,0])*np.power(10,iso_data[i,1])
        elif iso_data[i-1,0]>iso_data[i,0]:
            weights[i]=(iso_data[i+1,0]-iso_data[i,0])*np.power(10,iso_data[i,1])
        else:
            weights[i]=(iso_data[i+1,0]-iso_data[i-1,0])*np.power(10,iso_data[i,1])

    # Index data

    binned_data=[ [ [] for i in range(14) ] for i in range(14) ]

    for i in range(iso_data.shape[0]):

        try:
            binned_data[int(np.floor((iso_data[i,2]-3.3)/0.1))][int(np.floor((iso_data[i,3]+1.5)/0.5))].append(iso_data[:])
        except IndexError:
            print "Error:", iso_data[i,2], iso_data[i,3]



    # Re grid in hi-res

    metals=np.zeros(interp_points[0].shape)+metallicity

    
    logage=np.zeros(interp_points.T.shape[0])
    Mi=np.zeros(interp_points.T.shape[0])
    
    interp_photom={}
    for band in photom_data.keys():
        interp_photom[band] = np.zeros(interp_points.T.shape[0])
    
    for it, point in enumerate(interp_points.T):
        if it%10000==0:
#            print "\r", it, interp_points.shape[1],
            stdout.write("\r{0} of {1}".format(it,interp_points.shape[1]) )
            stdout.flush()

        selection_array=np.power(iso_data[:,2]-point[0],2)*4. + np.power(iso_data[:,3]-point[1],2)< 0.01

        shortlist_iso=iso_data[selection_array]
        shortlist_weights=weights[selection_array]
        shortlist_photom={}
        for band in photom_data.keys():
            shortlist_photom[band]=photom_data[band][selection_array]
        
        if shortlist_iso.size==0:
            logage[it]=np.nan
            Mi[it]=np.nan

            for band in bands_dict.values():
                interp_photom[band][it] = np.nan        

        else:
            short_weights1=shortlist_weights[:,-1]/(1E-9+np.power(shortlist_iso[:,2]-point[0],2.)*36.+np.power(shortlist_iso[:,3]-point[1],2.))
            short_weights_sum1=np.sum(short_weights1)
            short_weights2=1/(1E-9+np.power(shortlist_iso[:,2]-point[0],2.)*36.+np.power(shortlist_iso[:,3]-point[1],2.))
            short_weights_sum2=np.sum(short_weights2)

            logage[it]=np.sum(shortlist_iso[:,1]*short_weights1)/short_weights_sum1
            Mi[it]=np.sum(shortlist_iso[:,0]*short_weights1)/short_weights_sum1
            
            for band in bands_dict.values():
                interp_photom[band][it]=np.sum(shortlist_photom[band]*short_weights2)/short_weights_sum2
                
    inner_counts, xe, ye=np.histogram2d(iso_data[:,2], iso_data[:,3], bins=[logT_edges, logg_edges])
    outer_counts=(snd.filters.uniform_filter(inner_counts, size=(2,5))*25+0.1).astype(int)
    inner_Jac, xe, ye=np.histogram2d(iso_data[:,2], iso_data[:,3], bins=[logT_edges, logg_edges], weights=weights[:].flatten())
    outer_Jac=(snd.filters.uniform_filter(inner_Jac, size=(2,5))*25)

    outer_Jac[outer_Jac<0.01]=0.
    
    # Save to file
    
    output_array=[metals, Mi, logage, interp_points[0], interp_points[1], outer_Jac.flatten()]
    fmt_list=['%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3e']
    for band in bands_ordered:
        output_array.append(interp_photom[band])
        fmt_list.append('%.3f')
    output_array.extend([inner_counts.flatten(), outer_counts.flatten()])
    fmt_list.extend(['%i', '%i'])    

    np.savetxt(output_obj, np.array(output_array).T, fmt=fmt_list)
    
    
    
    
    
    
    
def padova_interpolated_isomake(directories, bands_dict, output_filename, bands_ordered=None):
    """ padova_interpolated_isomake(directories, bands_dict, 
                                    output_filename, 
                                    bands_ordered=None)

        Interpolate several isochrones defined in 
        (initial Mass, log(age)) space onto a grid of points
        in (log(T_eff), log(g)).

        This has been designed for use with Padova isochrones.

        Parameters
        ----------
        directories : string, list(string)
            the directories that contain the isochrone
            files
        bands_dict : dict
            A dictionary containing the photometric bands
            to be used.
        output_filename : 
            a filename object to save the interpolated 
            isochrones to.
        bands_ordered : list, optional
            The bands in their desired order. If not given
            is made from the values of bands_dict.
    """

    if isinstance(directories, basestring):
        directories=[directories]
        
    if bands_ordered is None:
        bands_ordered=bands_dict.values()        
        
    output_obj=open(output_filename, "w")
    
    header_string="#\t[M/H]\tMi\tlogAge\tlogTe\tlogg\tJacobian"
    for band in bands_ordered:
        header_string+="\t{}".format(band)
    header_string+="\tinner_count\touter_count\n"
    output_obj.write(header_string)
        
        

    iso_metal_dict={}        
    bands_metal_dicts={}
    for band in bands_dict.keys():
        bands_metal_dicts[band]={}
    
    # instead do this on band-by-band basis? *******************
        
    for direc in directories:
        iso_files_gz=gb.glob("{}/*.dat.gz".format(direc.rstrip("/")) )
        iso_files=gb.glob("{}/*.dat".format(direc.rstrip("/")) )
        
        # check for metallicity of each file
        # and check which bands it has

        for iso_file1 in iso_files_gz:
            metal=None
            iso_data=gz.open("{0}".format(iso_file1) )
            for line in iso_data:
                split_line=line.split()
                if "[M/H]" in split_line:
                    metal=float(split_line[split_line.index("[M/H]")+2])
                if "M_ini" in split_line:
                    for band in bands_metal_dicts.keys():
                        if band in split_line:
                            bands_metal_dicts[band][metal]=iso_file1
                                   
        for iso_file1 in iso_files:
            metal=None
            iso_data=open("{0}".format(iso_file1), "r" )
            for line in iso_data:
                split_line=line.split()
                if "[M/H]" in split_line:
                    metal=float(split_line[split_line.index("[M/H]")+2])
                if "M_ini" in split_line:
                    for band in bands_metal_dicts.keys():
                        if band in split_line:
                            bands_metal_dicts[band][metal]=iso_file1                      

        
    for metal in bands_metal_dicts[bands_metal_dicts.keys()[0]]:
        filenames=[ ]
        for band in bands_metal_dicts:
            if metal in bands_metal_dicts[band]:
                if bands_metal_dicts[band][metal] not in filenames:
                    filenames.append(bands_metal_dicts[band][metal])
            else:
                break
        else:
            iso_metal_dict[metal]=filenames
            
    print iso_metal_dict
    keys=iso_metal_dict.keys()
    keys.sort()

    if len(keys)>2:
#        iso_metal_weights=dict(zip(keys, np.gradient(np.array(keys)) ) )       in numpy 1.9.0 gradient has changed to use second order behaviour at boundaries which gives wrong results in this context
        iso_metal_weights=dict(zip(keys, replacement_gradient(np.array(keys)) ) ) 
    else:
        iso_metal_weights=dict(zip(keys, np.ones(len(keys)) ) )
    print "metals and weights: ",   iso_metal_weights      
    
# interp in metallicity order

    for key in keys:
	    iso_interp(iso_metal_dict[key], key, iso_metal_weights[key], output_obj, bands_dict, bands_ordered)

    output_obj.close()    
        
        
