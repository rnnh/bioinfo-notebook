# Using Ubuntu through a Virtual Machine

*Ubuntu* is a Linux operating system that is widely used for bioinformatics.
If you have not used a Linux system before, an Ubuntu virtual machine is an ideal way to try the programs documented on this website.

A *virtual machine* is a computer file, typically called an image, that behaves like an actual computer.
It acts as a computer within a computer.
Virtual machines run in a window, much like any other program running in a window on your computer.
The virtual machine is sequestered from the rest of the system, meaning that the software inside a virtual machine can not tamper with the computer itself.
This produces an ideal environment for testing other operating systems, and running software or applications on operating systems they were not originally intended for.

An Ubuntu virtual machine can be created using *VirtualBox*, and an Ubuntu *disk image*.
VirtualBox is a program that can be used to create, manage, and access virtual machines.
A disk image is a file that acts like a compact disc, or another storage device.
VirtualBox and the Ubuntu disk image are freely available online.

## Files required to set up an Ubuntu virtual machine

To set up an Ubuntu virtual machine, you will need an Ubuntu disk image, a file to install VirtualBox, and the VirtualBox Extension Package.
This requires around 13 GB of free hard drive space on your computer in total.
The Ubuntu disk image is around 2 GB in size, and may take a while to download depending on your internet connection.
The file required to install VirtualBox is around 108 or 123 MB in size, depending on the platform of your computer (i.e. Windows or Mac).

### Direct links to download required files

1. [The Ubuntu disk image (filename: `ubuntu-18.04.4-desktop-amd64.iso`)](http://releases.ubuntu.com/18.04.4/ubuntu-18.04.4-desktop-amd64.iso)
2. [VirtualBox installer for Windows](https://download.virtualbox.org/virtualbox/6.1.4/VirtualBox-6.1.4-136177-Win.exe)
3. [VirtualBox installer for Mac](https://download.virtualbox.org/virtualbox/6.1.4/VirtualBox-6.1.4-136177-OSX.dmg)
4. [VirtualBox Extension Pack (all platforms)](https://download.virtualbox.org/virtualbox/6.1.4/Oracle_VM_VirtualBox_Extension_Pack-6.1.4.vbox-extpack)

If the above links do not work, they may have expired.
In this case, the above files can be found on the [VirtualBox website](https://www.virtualbox.org/wiki/Downloads) and the [Ubuntu website](https://ubuntu.com/download/desktop).

## How to create an Ubuntu virtual machine using VirtualBox

1. Download the [VirtualBox installer](#direct-links-to-download-required-files) for your computer (either Windows or Mac).
2. Once the VirtualBox installer is downloaded, open it and follow the on-screen instructions to install the VirtualBox program.
3. **Windows only:** If you get a "Windows Security" prompt asking *"Would you like to install this device software?"* for driver software from *"Publisher: Oracle Corporation"*, select "Install".
4. **Mac only:** If you get a *"This package will run a program to determine if the software can be installed"* prompt while installing VirtualBox, select "Continue". You may also be asked to enter your user password while installing VirtualBox on a Mac.
5. Once installed, open the VirtualBox program.
6. In VirtualBox, click on "New" (the blue badge). This will open a menu to create a new virtual machine.
7. In the "Name" field of the "Name and operating system" window, type "ubuntu". VirtualBox will automatically set the type and version for this virtual machine as "Linux" and "Ubuntu".
8. Select "Next" to proceed to the "Memory size" section.
9. In this section, you can set the amount of Random Access Memory (RAM) that the virtual machine can use. A suggested amount of RAM will automatically be selected when you get to this page, but you can increase the amount of RAM allocated using the slider on this page.
10. **Note:** If you use the slider to increase the amount of RAM allocated on the "Memory Size" page, keep the slider in the green zone. Setting the slider in the orange or red zone (>50% of your computer's available RAM) will negatively affect the performance of the virtual machine.
11. Select "Next" to proceed to the "Hard disk" page.
12. Select "Create a virtual hard disk now", and then select "Create".
13. On the "Hard disk file type" page, select "VDI (VirtualBox Disk Image)", and then select "Next".
14. On the "Store on physical hard disk", select "Dynamically allocated", and then select "Next" to proceed to the "File location and size" page.
15. On this page, you can change the location and size of the virtual hard disk. There is no need to adjust the size of the virtual hard disk, but take note of its location (the folder/directory it will be created in). Select "Create".
16. In the left side of the VirtualBox main menu, double-click the name of the virtual machine you just created ("ubuntu").
17. This will bring up the "Select start-up disk" window. In this window, select the folder icon to open the "Optical Disk Selector" menu.
18. In this menu, select "Add", which will open a window titled "Please choose a virtual optical disk file".
19. In this window, go to the folder into which the Ubuntu disk image downloaded (e.g. "Downloads"), and click the [Ubuntu disk image (filename: `ubuntu-18.04.4-desktop-amd64.iso`)](#direct-links-to-download-required-files) to select it, and then select "Open".
20. This will bring you back to the "Optical Disk Selector" window. Select the Ubuntu disk image you selected in the previous window, and click on "Choose".
21. This will bring you back to the "Select start-up disk" window. The Ubuntu disk image should be selected in the drop down menu (this read "Empty" before the Ubuntu disk image was added). Select "Start" to start the virtual machine.
22. The Ubuntu virtual machine is now running in its own window. It may take a few minutes to start up the first time.
23. On the "Welcome" screen in Ubuntu, select "Install Ubuntu".
24. In the "Keyboard layout" section, select your keyboard layout, and then select "Continue". This will bring you to the "Updates and other software" window.
25. In this window, in the section "What apps would you like to install to start with?", select "Minimal installation".
26. In the "Other options" section, select "Download updates while installing Ubuntu", and leave "Install third-party software..." unselected. Select "Continue" to proceed to the "Installation type" window.
27. In this window, select "Erase disk and install Ubuntu". As this is a virtual machine, in this instance "disk" refers to the virtual disk image (`.vdi`) file created earlier (see steps 12 to 15). Select "Install now".
28. A window titled "Write the changes to disks?" will appear. In this window, select "Continue".
29. This will bring you to the "Where are you?" window. In this window, enter your location (which is needed to set the system clock) and select "Continue".
30. Fill in the requested details in the "Who are you?" window: your name, your computer's name, your username (both of which will be filled in automatically when you enter your name), and your password. Make sure you remember your password, you will need it to install programs in your Ubuntu virtual machine. Select "Continue" to proceed.
31. At this point, Ubuntu will begin installing on the virtual disk image created earlier (the `.vdi` file). This will take a few minutes.
32. Once the installation is complete, select "Restart Now" from the "Installation complete" dialog window.
33. When asked "Please remove the installation media and press ENTER", press Enter (a.k.a. Return).
34. The virtual machine will then restart, and the Ubuntu login page will load. On this page, select the user you created during the installation, and enter your password to log in.
35. Once you have logged in, you have finished setting up your Ubuntu virtual machine. Click through the "What's new in Ubuntu" window for a brief introduction to Ubuntu.
36. When you want to close your Ubuntu virtual machine, close the window it is running in to bring up the "Close Virtual Machine" window, select "Power off the machine" and click "OK". This is the equivalent of shutting down the machine. Alternatively, you can select "Power off" within the Ubuntu virtual machine.

Once you have finished installing the Ubuntu virtual machine, you can delete the Ubuntu disk image (filename: `ubuntu-18.04.4-desktop-amd64.iso`), and the VirtualBox installer.

## Increasing the screen resolution of the Ubuntu virtual machine

At this point, the Ubuntu virtual machine takes up only a small portion of the VirtualBox window it runs in.
To increase the screen resolution of the Ubuntu virtual machine, you will need to download the [VirtualBox Extension Package](#direct-links-to-download-required-files) and follow the steps below.

1. Once downloaded, double click the VirtualBox Extension Pack (file extension `.vbox-extpack`). If you have installed the VirtualBox program, it will open this file.
2. VirtualBox will open with a window notifying that an extension pack is about to be installed. In this window, select "Install" to proceed with the extension pack installation.
3. Scroll to the bottom of the Terms and Conditions window that opens, and select "I Agree" to install the extension pack.
4. Open the Ubuntu virtual machine in VirtualBox.
5. In the menu bar of the VirtualBox window in which Ubuntu is running, select the "Devices" menu, and select "Insert Guest Additions CD image...".
6. A notification will appear in the Ubuntu virtual machine: '"VBox_GAs_6.1.4" contains software intended to be automatically started. Would you like to run it?". In this window, select "Run", and enter your Ubuntu password to install the VirtualBox Guest Additions on the Ubuntu virtual machine.
7. A terminal window will open showing the VirtualBox Guest Additions installation progress. Once the installation has finished, press Return (Enter) to close this window.
8. Close the Ubuntu virtual machine by closing the window it is running in, and selecting "Power off the machine" from the "Close Virtual Machine" window.
9. Open the Ubuntu virtual machine in VirtualBox.
10. In the menu bar of the window in which Ubuntu is running, select the "View" menu, and confirm that "Auto-resize Guest Display" is enabled.

## See also

- [Windows Subsystem for Linux](wsl.md)
- [The Ubuntu Website](https://ubuntu.com/)
- [The VirtualBox Website](https://www.virtualbox.org/)

## References

- [What is a Virtual Machine?](https://azure.microsoft.com/en-us/overview/what-is-a-virtual-machine/)
- [How to Install Ubuntu on VirtualBox](https://www.wikihow.com/Install-Ubuntu-on-VirtualBox)
