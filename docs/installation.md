---
title: Installation and Dependencies
sidebar_position: 3
---

## Installation from the GitHub repository

The easiest way to install TidyScreen is to download the installation script from the [repository](https://github.com/alfredoq/TidyScreen_v2). 
You can start working after downloading the bash shell script and executing as follows:

```bash
chmod 777 tidyscreen_installation.sh
./tidyscreen_installation.sh

## The creation con dedicated environment named 'tidyscreen' will be accomplished

conda activate tidyscreen

# Start working
```

Upon execution of the bash shell script, two new `conda`environments will be created:
- `tidyscreen`
- `adt`

In order to work, you only need to activate the `tidyscreen` environment, since the one corresponding `adt` will only be used by TidyScreen internally for specific operations.


## Manual installation of the package and its dependencies

Follow the instructions provided en README of the [repository](https://github.com/alfredoq/TidyScreen_v2).