import mrcfile
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler

# File path to the MRC file
ptcls_file = "/home/elmlundho/data/NanoX/dissolving/NP87/subtr.mrc"
# Save the denoised particle stack to a new MRC file in Downloads/nano
file_path1  = "/home/elmlundho/data/NanoX/dissolving/NP87/subtr_den.mrc"  # Change 'your-username' to your actual username

# Open and update the header of the MRC file
with mrcfile.open(ptcls_file, mode='r+', permissive=True) as mrc:
    mrc.update_header_from_data()

# Read the particle stack from the MRC file
ptcls_stk = mrcfile.read(ptcls_file)

# Get the shape of the particle stack
stk_shape = np.shape(ptcls_stk)
ptcls_no = stk_shape[0]  # Number of particle images in the stack

# No batching
num_batches = 1

dnmat = np.zeros((ptcls_no, stk_shape[1], stk_shape[2]))

# First pass of denoising with RBF kernel
for j in range(num_batches):
    start_index = 0
    end_index =  ptcls_no
    current_batch_size = end_index - start_index
    mat = np.zeros((current_batch_size, stk_shape[1], stk_shape[2]))
    
    for i in range(current_batch_size):
        mat[i, :, :] = ptcls_stk[start_index + i, :, :]

    flattened_images = mat.reshape(current_batch_size, -1)

    # Optionally standardize the data to have mean zero and variance one
    scaler = StandardScaler()
    flattened_images_std = scaler.fit_transform(flattened_images)

    # Apply KernelPCA
    kpca = KernelPCA(n_components=500, kernel="cosine", fit_inverse_transform=True)
    transformed_data = kpca.fit_transform(flattened_images_std)
    denoised_data = kpca.inverse_transform(transformed_data)
    denoised_data_rescaled = scaler.inverse_transform(denoised_data)

    # Reshape the denoised data back to the original image shape
    dn = denoised_data_rescaled.reshape(current_batch_size, stk_shape[1], stk_shape[2])

    for i in range(current_batch_size):
        dnmat[start_index + i, :, :] = dn[i, :, :]
    
dnf = dnmat.astype(np.float32)
with mrcfile.new(file_path1, overwrite=True) as mrc:	
    mrc.set_data(dnf)



