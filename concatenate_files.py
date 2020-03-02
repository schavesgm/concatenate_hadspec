"""
    Written by Sergio Chaves Garcia-Mascaraque as a modification of
    a code by Benjamin Jaegger.
"""

import numpy as np
import pandas as pd
import sys, os, re

def read_directory(directory):
    """
        Read the files in a directory assumed to collect openQCD hadspec files
        and output a pandas dataframe with the data

        Parameters
        ----------
        directory : string
            name of directory where the files are located

        Returns
        -------
        data : pandas.DataFrame
            Dataframe with the following columns:
                - file : full filename (without directory)
                - run : the openQCD run name
                - config : the configuration number
                - source : the source number
                - particle : a string particle identifier
                - flavour : the flavour combination of the particle
    """

    column_names = ["file", "run", "config", "source", "particle", "flavour"]
    filelist = []

    for file in os.listdir(directory):

        m = re.match(r'(.*?)n(\d+)\.s(\d+)\.(.*)\.(.*)', file)

        if not m:
            continue

        filelist.append([
            str(m.group(0)),
            str(m.group(1)),
            int(m.group(2)),
            int(m.group(3)),
            str(m.group(4)),
            str(m.group(5))
        ])


    data = pd.DataFrame( filelist, columns=column_names )
    data["run"] = data["run"].astype("category")
    data["particle"] = data["particle"].astype("category")
    data["flavour"] = data["flavour"].astype("category")

    return data

def print_meson_file( sub_dataFrame, in_dir, out_dir, Ng = 16 ):
    """
        Concat the contents of meson files and store them to disk, there will be
        16 output files, one for each gamma matrix

        Parameters
        ----------
        subDataFrame : pandas.DataFrame
            dataframe produced by read_directory filtered to be for one particle
            and flavour combination

        in_dir : string
            The directory the files to be read are in

        out_dir : string
            The directory to store the concat results

        Ng      : Integer ( optional = 16 )
            Total number of gamma structures calculated.
    """

    # Get the first row of the data frame
    row1 = sub_dataFrame.iloc[0]

    # Load the data corresponding to get the NT of the file
    file_values = np.loadtxt( "{}/{}".format( in_dir, row1["file"]) )
    nt = file_values.shape[0]
    numPsquare = int( file_values.shape[1] / ( 2 * Ng ) )

    # Create a folder to hold the data in an organised way
    folder_partflav = "{}/{}".format( row1["particle"], row1["flavour"] )

    # Generate string formats to be filled in the generation of data
    format_file = "{}.{}.g{}.{}.p{}".format(
                row1["run"], row1["particle"], "{}", row1["flavour"], "{}" )
    format_gamma = "{}/{}/g_{}".format( out_dir, folder_partflav, "{}" )
    format_momenta = "{}/psquare_{}".format( "{}", "{}" )

    # Buffer to save the data to, consisting of nt + 16 * ( Re + Im )
    buff_data = np.zeros( ( nt, 1 + numPsquare * 2 * Ng ) )
    buff_data[:, 0] = range(nt)

    files_to_save = []
    for ig in range( Ng ):
        path_gamma = format_gamma.format( ig )
        if not os.path.exists( path_gamma ):
            os.makedirs( path_gamma )
        for ip in range( numPsquare ):
            path_momenta = format_momenta.format( path_gamma, ip )
            if not os.path.exists( path_momenta ):
                os.makedirs( path_momenta )
            # Create the file inside the gamma/momenta directory
            path_file = "{}/{}".format( path_momenta, format_file.format( ig, ip ) )
            files_to_save.append( open( path_file, "w" ) )

    for index, row_db in sub_dataFrame.iterrows():

        # Load the data corresponding to each config and source into the buffer
        buff_data[:, 1:] = np.loadtxt( "{}/{}".format( in_dir, row_db["file"] ) )

        # Append to the correspoding file with coordinates ( ig, ip )
        for ig in range( Ng ):
            for ip in range( numPsquare ):
                index_real = 2 * ( ig * numPsquare + ip ) + 1
                index_imag = 2 * ( ig * numPsquare + ip ) + 2
                np.savetxt( files_to_save[ig*numPsquare + ip],
                            buff_data[:, [0, index_real, index_imag] ],
                            fmt="%d\t%.12e\t%.12e" )

    for files in files_to_save:
        files.close()

def print_baryon_file( sub_dataFrame, in_dir, out_dir ):
    """
        Concat the contents of meson files and store them to disk, there will be
        16 output files, one for each gamma matrix

        Parameters
        ----------
        subDataFrame : pandas.DataFrame
            dataframe produced by read_directory filtered to be for one particle
            and flavour combination

        in_dir : string
            The directory the files to be read are in

        out_dir : string
            The directory to store the concat results

        Ng      : Integer ( optional = 16 )
            Total number of gamma structures calculated.
    """

    # Get the first row of the data frame
    row1 = sub_dataFrame.iloc[0]

    # Load the data corresponding to get the NT of the file
    file_values = np.loadtxt( "{}/{}".format( in_dir, row1["file"]) )
    nt = file_values.shape[0]

    # Number of momenta in the calculation
    numPsquare = int( file_values.shape[1] / ( 2 ) )

    # Create a folder to hold the data in an organised way
    folder_partflav = "{}/{}".format( row1["particle"], row1["flavour"] )

    # Generate string formats to be filled in the generation of data
    format_file = "{}.{}.{}.{}.p{}".format(
            row1["run"], row1["particle"], row1["particle"], row1["flavour"], "{}" )
    format_momenta = "{}/baryon/{}/psquare_{}".format( out_dir, folder_partflav, "{}" )

    # Buffer to save the data to, consisting of nt + 16 * ( Re + Im )
    buff_data = np.zeros( ( nt, 1 + numPsquare * 2  ) )
    buff_data[:, 0] = range(nt)
    
    # Get all the files that are going to be analysed
    files_to_save = []
    for ip in range( numPsquare ):
        # Path to a given momenta
        path_momenta = format_momenta.format( ip )
        if not os.path.exists( path_momenta ):
            os.makedirs( path_momenta )
        # Create the file inside the momenta directory
        path_file = "{}/{}".format( path_momenta, format_file.format( ip ) )
        files_to_save.append( open( path_file, "w" ) )

    for index, row_db in sub_dataFrame.iterrows():

        # Load the data corresponding to each config and source into buffer
        buff_data[:, 1:] = np.loadtxt( "{}/{}".format( in_dir, row_db["file"] ) )

        # Append to the corresponding file with coordinate ip
        for ip in range( numPsquare ):
            index_real = 2 * ip + 1
            index_imag = 2 * ip + 2
            np.savetxt( files_to_save[ip],
                        buff_data[:, [0, index_real, index_imag] ],
                        fmt="%d\t%.12e\t%.12e" )

    for files in files_to_save:
        files.close()

def main(args):
    """
        Run the script to concatenate the files for each different flavour and
        particle.

        Parameters
        ----------
        args[1] : string
            The input directory, location of the files

        args[1] : string ( optional = args[1] )
            The ouput directory where the concatenated files will be located
    """

    in_dir = args[1]

    if len( args ) < 3:
        out_dir = in_dir
    else:
        out_dir = args[2]

    if len( args ) == 4:
        Ng = args[3]
    else:
        Ng = 16
    
    # Read the data inside the input directory
    dataFrame = read_directory( in_dir )

    if len( dataFrame["run"].cat.categories ) > 1:
        print( "Sorting multiple runs currently not implemented" )
        return 1

    # Loop over different flavour combinations
    for flavour in dataFrame["flavour"].cat.categories:
        flavor_dataFrame = dataFrame[ dataFrame["flavour"] == flavour ]

        # Loop over different particles
        for particle in flavor_dataFrame["particle"].cat.categories:
            particle_dataFrame = flavor_dataFrame[
                                flavor_dataFrame["particle"] == particle]

            # Skip empty selections
            if ( len( particle_dataFrame ) ) == 0:
                continue

            if particle == "meson":
                print_meson_file( particle_dataFrame, in_dir, out_dir, Ng )
            else:
                print_baryon_file( particle_dataFrame, in_dir, out_dir )

if __name__ == '__main__':
    try:
        main( sys.argv )
    except KeyboardInterrupt:
        print( "\nProgram interrupted by user" )
