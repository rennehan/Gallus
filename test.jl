import Gallus

using DataFrames, Printf

#num_halos, num_galaxies = read_file(FileTypes.Rockstar("halos_0.0.ascii"))

#@printf("Found %d halos and %d galaxies.\n\n", num_halos, num_galaxies)

#@printf("The first halo ID is %d\n", Halos.master_list[1].halo_id)
#@printf("The first halo mass is %g Msun/h\n\n", Halos.master_list[1].masses["mvir"])

data_directory = "/Users/rennehan/Downloads"
s = Simulation(FileTypes.GIZMO("$data_directory/snapshot_169.0.hdf5", "$data_directory/MFM_Simba_Vida_N128L12.tex"))
