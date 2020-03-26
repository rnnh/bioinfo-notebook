# Using Ubuntu throught a Virtual Machine

## How to create an Ubuntu virtual machine using VirtualBox

Note, this requires ~13 GB of free hard drive space. The `ubuntu-18.04.4-desktop-amd64.iso` disk image is ~2 GB. The installers for VirtualBox require 108 to 123 MB.

1. To emulate Ubuntu, a disk image (also known as an ISO file, with the file extension `.iso`) is needed. An Ubuntu disk image can be downloaded using this link: <http://releases.ubuntu.com/18.04.4/ubuntu-18.04.4-desktop-amd64.iso>
2. VirtualBox is also needed to emulate Ubuntu. Follow this link for the VirtualBox downloads page: <https://www.virtualbox.org/wiki/Downloads>
3. On the VirtualBox downloads page, select the platform package for your computer (in the "platform packages" section). If you are using Windows, select "Windows hosts"; if you're using a Mac, select "OS X hosts".
4. Once the VirtualBox platform package is downloaded, open it and follow the on-screen instructions to install the VirtualBox program.
5. **Windows only:** If you get a "Windows Security" prompt asking *"Would you like to install this device software?"* for driver software from *"Publisher: Oracle Corporation"*, select "Install".
6. **Mac only:** If you get a *"This package will run a program to determine if the software can be installed"* prompt while installing VirtualBox, select "Continue". You may also be asked to enter your user password while installing VirtualBox on a Mac.
7. Once installed, open the VirtualBox program.
8. In VirtualBox, click on "New" (the blue badge). This will open a menu to create a new virtual machine.
9. In the "Name" field of the "Name and operating system" window, type "ubuntu". VirtualBox will automatically set the type and version for this virtual machine as "Linux" and "Ubuntu".
9. Select "Next" to proceed to the "Memory size" section.
10. In this section, you can set the amount of Random Access Memory (RAM) that the virtual machine can use. A suggested amount of RAM will automatically be selected when you get to this page, but you can increase the amount of RAM allocated using the slider on this page.
11. **Note:** If you use the slider to increase the amount of RAM allocated on the "Memory Size" page, keep the slider in the green zone. Setting the slider in the orange or red zone will negatively affect the performance of the virtual machine.
12. Select "Next" to proceed to the "Hard disk" page.
13. Select "Create a virtual hard disk now", and then select "Create".
14. On the "Hard disk file type" page, select "VDI (VirtualBox Disk Image)", and then select "Next".
15. On the "Store on physical hard disk", select "Dynamically allocated", and then select "Next" to proceed to the "File location and size" page.
16. On this page, you can change the location and size of the virtual hard disk. There is no need to adjust the size of the virtual hard disk, but take note of its location (the folder/directory it will be created in). Select "Create".
17. In the left side of the VirtualBox main menu, double-click the name of the virtual machine you just created ("ubuntu").
18. This will bring up the "Select start-up disk" window. In this window, select the folder icon to open the "Optical Disk Selector" menu.
19. In this menu, select "Add", which will open a window titled "Please choose a virtual optical disk file".
20. In this window, go to the folder into which the Ubuntu disk image downloaded (e.g. "Downloads"), and click the Ubuntu disk image file (ending in `.iso`) to select it, and then select "Open".
21. This will bring you back to the "Optical Disk Selector" window. Select the Ubuntu disk image you selected in the previous window, and click on "Choose".
22. This will bring you back to the "Select start-up disk" window. The Ubuntu disk image should be selected in the drop down menu (this read "Empty" before the Ubuntu disk image was added). Select "Start" to start the virtual machine.
23. The Ubuntu virtual machine is now running in its own window. It may take a few minutes to start up the first time.
24. On the "Welcome" screen in Ubuntu, select "Install Ubuntu".
25. In the "Keyboard layout" section, select your keyboard layout, and then select "Continue". This will bring you to the "Updates and other software" window.
26. In this window, in the section "What apps would you like to install to start with?", select "Minimal installation".
27. In the "Other options" section, select "Download updates while installing Ubuntu", and leave "Install third-party software..." unselected. Select "Continue" to proceed to the "Installation type" window.
28. In this window, select "Erase disk and install Ubuntu". As this is a virtual machine, in this instance "disk" refers to the virtual disk `.vdi` file created earlier (see steps 13 to 16). Select "Install now".
29. A window titled "Write the changes to disks?" will appear. In this window, select "Continue".
30. This will bring you to the "Where are you?" window. In this window, enter your location (which is needed to set the system clock) and select "Continue".
31. Fill in the requested details in the "Who are you?" window: your name, your computer's name, your username (both of which will be filled in automatically when you enter your name), and your password. Make sure you remember your password, you will need it to install programs in your Ubuntu virtual machine. Select "Continue" to proceed.
32. At this point, Ubuntu will begin installing on the virtual disk image created earlier (the `.vdi` file). This will take a few minutes.
33. Once the installation is complete, select "Restart Now" from the "Installation complete" dialog window.
34. When asked "Please remove the installation media and press Enter", press Enter.
35. The virtual machine will then restart, and the Ubuntu login page will load. On this page, select the user you created during the installation, and enter your password to log in.
36. Once you have logged in, you have finished setting up your Ubuntu virtual machine. Click through the "What's new in Ubuntu" window for a brief introduction to Ubuntu.
37. When you want to close your Ubuntu virtual machine, close the window it is running in to bring up the "Close Virtual Machine" window, select "Power off the machine" and click "OK". This is the equivalent of shutting down the machine. Alternatively, you can select "Power off" within the Ubuntu virtual machine.

Once you have finished installing the Ubuntu virtual machine, you can delete the disk image downloaded in step 1.

## Increasing the screen resolution of the Ubuntu virtual machine

1. In the "VirtualBox Manager" window, select the "Ubuntu" virtual machine, and select "Settings" (the gear icon next to "New").
2. In the left hand menu of "Settings", select "Display".
3. In the "Screen" tab of "Display", increase "Video memory" to its maximum by moving the slider to the right.
4. In the "Screen" tab of "Display", select "Enable 3D Acceleration" and "Enable 2D Acceleration", where available.
5. Click on "OK" to apply these settings.
6. Start the Ubuntu virtual machine.
7. In the menu bar of the virtual box window, select the "Devices" menu, and select "Insert Guest Additions CD image...".
8. On the desktop in the Ubuntu virtual machine, double click the "VBox_GAs_..." disk icon.
9. This will open the disk in a window; select "Run Software" from the dialog at the top of the window.
10. Select "Run" from the dialog box that appears, and enter your Ubuntu password.
11. After restarting the Ubuntu virtual machine, it will now scale to the full size of the window it is in.