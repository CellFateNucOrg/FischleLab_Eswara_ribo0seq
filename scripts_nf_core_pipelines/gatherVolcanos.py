import os
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Set the root directory
root_dir = '/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq'
sub_dir = 'plots/differential'
samples = os.listdir(os.path.join(root_dir, sub_dir))

images_info = []

for sample in samples:
    image_dir = os.path.join(root_dir, sub_dir, sample,'png')
    if(os.path.isdir(os.path.join(image_dir))):
        png_file = [f for f in os.listdir(image_dir) if f.lower().endswith('.png')]
        images_info.append((sample, os.path.join(root_dir,sub_dir,sample,'png',png_file[0])))


# Paginate in chunks of 6 (3x2)
chunks = [images_info[i:i+6] for i in range(0, len(images_info), 6)]

# Create multipage PDF
with PdfPages(os.path.join(root_dir,'custom/plots/volcanos.pdf')) as pdf:
    for chunk in chunks:
        fig, axs = plt.subplots(2, 3, figsize=(11, 8.5))
        axs = axs.flatten()
        for ax, (folder, img_path) in zip(axs, chunk):
            img = Image.open(img_path)
            ax.imshow(img)
            ax.set_title(os.path.basename(folder), fontsize=10)
            ax.axis('off')
        for ax in axs[len(chunk):]:
            ax.axis('off')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
