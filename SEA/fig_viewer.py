import os
import re
import matplotlib.pyplot as plt
from ipywidgets import interact, IntSlider

SEA_path = os.path.dirname(os.path.abspath(__file__))


exp_path = "exp_2023-10-29-12-10-02"
figs_legvel_path = f"{SEA_path}/{exp_path}/figs_legvel"

# Get the list of image files
image_files = sorted(os.listdir(figs_legvel_path))

# Extract the values of m, t, and L from the image file names
# m_values = sorted(list(set(re.findall(r"_m(\d+)_", " ".join(image_files)))))
# t_values = sorted(list(set(re.findall(r"_t(\d+)_", " ".join(image_files)))))
# L_values = sorted(list(set(re.findall(r"_L(\d+)", " ".join(image_files)))))


m_values = sorted([int(m) for m in set(re.findall(r"_m(\d+)_", " ".join(image_files)))])
t_values = sorted([int(t) for t in set(re.findall(r"_t(\d+)_", " ".join(image_files)))])
L_values = sorted([int(L) for L in set(re.findall(r"_L(\d+)", " ".join(image_files)))])


# m_values = sorted([int(m) for m in sorted(list(set(re.findall(r"_m(\d+)_", " ".join(image_files)))))])
# t_values = sorted([int(t) for t in sorted(list(set(re.findall(r"_t(\d+)_", " ".join(image_files)))))])
# L_values = sorted([int(L) for L in sorted(list(set(re.findall(r"_L(\d+)", " ".join(image_files)))))])






print(L_values)

@interact(m=IntSlider(min=min(m_values), max=max(m_values), step=1, value=min(m_values)),
          t=IntSlider(min=min(t_values), max=max(t_values), step=1, value=min(t_values)),
          L=IntSlider(min=min(L_values), max=max(L_values), step=1, value=min(L_values)))
def display_image(m, t, L):
    # Get the corresponding image file name
    image_file = f"k1_m{m}_t{t}_L{L}.png"

    # Check if the image file exists
    if image_file in image_files:
        # Load and display the image
        image_path = os.path.join(figs_legvel_path, image_file)
        image = plt.imread(image_path)
        plt.imshow(image)
        plt.axis('off')
        plt.show()
    else:
        print("Image not found.")