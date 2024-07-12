import mrcfile
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler

# File path to the MRC file
ptcls_file = "/Users/Polo/Downloads/nano/raw/NP87_raw.mrc"

# Open and update the header of the MRC file
with mrcfile.open(ptcls_file, mode='r+', permissive=True) as mrc:
    mrc.update_header_from_data()

# Read the particle stack from the MRC file
ptcls_stk = mrcfile.read(ptcls_file)

# Get the shape of the particle stack
stk_shape = np.shape(ptcls_stk)
ptcls_no = stk_shape[0]  # Number of particle images in the stack

# Determine the batch size to keep it around 90-110 particles
min_batch_size = 90
max_batch_size = 110
avg_batch_size = (min_batch_size + max_batch_size) // 2
num_batches = ptcls_no // avg_batch_size
remainder = ptcls_no % avg_batch_size

if remainder > 0:
    num_batches += 1
    avg_batch_size = ptcls_no // num_batches
    remainder = ptcls_no % num_batches

dnmat = np.zeros((ptcls_no, stk_shape[1], stk_shape[2]))

# First pass of denoising with RBF kernel
for j in range(num_batches):
    start_index = j * avg_batch_size
    end_index = start_index + avg_batch_size if j < num_batches - 1 else ptcls_no
    current_batch_size = end_index - start_index
    mat = np.zeros((current_batch_size, stk_shape[1], stk_shape[2]))
    
    for i in range(current_batch_size):
        mat[i, :, :] = ptcls_stk[start_index + i, :, :]

    flattened_images = mat.reshape(current_batch_size, -1)

    # Optionally standardize the data to have mean zero and variance one
    scaler = StandardScaler()
    flattened_images_std = scaler.fit_transform(flattened_images)

    # Apply KernelPCA
    kpca = KernelPCA(n_components=4, kernel="rbf", fit_inverse_transform=True)
    transformed_data = kpca.fit_transform(flattened_images_std)
    denoised_data = kpca.inverse_transform(transformed_data)
    denoised_data_rescaled = scaler.inverse_transform(denoised_data)

    # Reshape the denoised data back to the original image shape
    dn = denoised_data_rescaled.reshape(current_batch_size, stk_shape[1], stk_shape[2])

    for i in range(current_batch_size):
        dnmat[start_index + i, :, :] = dn[i, :, :]

# Second pass of denoising with cosine kernel
dnmat2 = np.zeros((ptcls_no, stk_shape[1], stk_shape[2]))
for j in range(num_batches):
    start_index = j * avg_batch_size
    end_index = start_index + avg_batch_size if j < num_batches - 1 else ptcls_no
    current_batch_size = end_index - start_index
    mat = np.zeros((current_batch_size, stk_shape[1], stk_shape[2]))
    
    for i in range(current_batch_size):
        mo = (start_index + i + avg_batch_size // 2) % ptcls_no
        mat[i, :, :] = dnmat[mo, :, :]

    flattened_images = mat.reshape(current_batch_size, -1)

    # Optionally standardize the data to have mean zero and variance one
    scaler = StandardScaler()
    flattened_images_std = scaler.fit_transform(flattened_images)

    # Apply KernelPCA with cosine kernel
    kpca = KernelPCA(n_components=4, kernel="cosine", fit_inverse_transform=True)
    transformed_data = kpca.fit_transform(flattened_images_std)
    denoised_data = kpca.inverse_transform(transformed_data)
    denoised_data_rescaled = scaler.inverse_transform(denoised_data)

    # Reshape the denoised data back to the original image shape
    dn = denoised_data_rescaled.reshape(current_batch_size, stk_shape[1], stk_shape[2])

    for i in range(current_batch_size):
        dnmat2[start_index + i, :, :] = dn[i, :, :]

# Save the denoised particle stack to a new MRC file in Downloads/nano
file_path = "/Users/Polo/Downloads/nano/NP87_purekPCAcosDN_NP4.mrc"  # Change 'your-username' to your actual username
dnf = dnmat2.astype(np.float32)
with mrcfile.new(file_path, overwrite=True) as mrc:
    mrc.set_data(dnf)
