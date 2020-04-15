---
layout: default
title: Windows Subsystem for Linux
parent: 1. General guides
nav_order: 2
---

# Windows Subsystem for Linux

Windows Subsystem for Linux (WSL) is a feature of Windows 10.
When enabled, WSL allows Linux systems (e.g. Ubuntu) to be used as Windows applications.
These Linux systems can be downloaded directly from the Microsoft Store.
The bioinfo-notebook [conda](conda.md) environment can be installed in an Ubuntu system running using WSL.

## Installing Ubuntu on Windows 10 using WSL

Before you begin, make sure you have around 1.20 GB of free disk space.

### Enable WSL

1. In the search box on the taskbar, type "control panel", and then select Control Panel.
2. In the Control Panel, select "Programs".
3. Under Programs and Features, select "Turn Windows features on or off".
4. If asked "Do you want this app to make changes to your device?", select "Yes".
5. From the list of Windows features, tick the box next to "Windows Subsystem for Linux" to enable WSL, and click OK.

### Download Ubuntu from the Microsoft Store

1. In the search box on the taskbar, type "microsoft store", and select Microsoft Store.
2. In the Microsoft Store, search for "Ubuntu".
3. Select the Ubuntu app.
4. On the app page, select "Get" to download Ubuntu.
5. If asked to sign in with a Microsoft account, select "No, thanks".

After enabling WSL and downloading Ubuntu from the Microsoft Store, Ubuntu can be used like a regular Windows application.

### Running Ubuntu for the first time

1. In the search box on the taskbar, type "Ubuntu", and select the Ubuntu app to launch it. It will take a few minutes to install the first time it runs.
2. When prompted, enter a UNIX username- this does not need to be the same as your Windows account name.
3. You will need to set a UNIX password. This is only used for the Ubuntu app, it does not need to be the same as your Windows password. Make sure you remember your UNIX password, as you will need it for installing new programs in Ubuntu.

Once your UNIX password has been updated successfully, you will see the `bash` command prompt in the Ubuntu window:

```
(Your UNIX username)@(Your computer's alias):~$ _
```

In this command prompt, the tilde character (`~`) indicates that you are currently in your home directory.
The dollar sign (`$`) indicates that this command line uses the `bash` shell language.

## See also

- [Introduction to the command line](cl_intro.md)
- [Using Ubuntu through a Virtual Machine](ubuntu_virtualbox.md) 
- [conda](conda.md)
