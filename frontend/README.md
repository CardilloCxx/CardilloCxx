# Install the required version
sudo apt update
sudo apt install gcc-13 g++-13

# Set as default using update-alternatives
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100

# Verify installation
gcc --version
g++ --version

# Install Vulkan Devtools
sudo apt install libvulkan-dev vulkan-tools

# Required Ubuntu dependency
sudo apt install xorg-dev libglu1-mesa-dev