#include "hdf5.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>

int main(void) {

  hid_t       file_id, dataset_id,dataspace_id; /* identifiers */
  herr_t      status;
  int         i, j, dset_data[4][6], read_data[4][6];
  hsize_t     dims[2];

  /* Initialize the dataset. */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 6; j++)
      dset_data[i][j] = i * 6 + j + 1;

  /* Create a new file using default properties. */
  file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Create the data space for the dataset. */
  dims[0] = 4;
  dims[1] = 6;
  dataspace_id = H5Screate_simple(2, dims, NULL);

  /* Create the dataset. */
  dataset_id = H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write the dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    dset_data);

  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);

//------------------------------------------------------

  /* Open an existing dataset. */
  dataset_id = H5Dopen2(file_id, "/dset", H5P_DEFAULT);


  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   read_data);

  for (i = 0; i < 4; i++)
    for (j = 0; j < 6; j++)
      printf("%d ",read_data[i][j]); // 1-24

  /* Close the dataset. */
  status = H5Dclose(dataset_id);

  /* Close the file. */
  status = H5Fclose(file_id);

  return 0;
}
