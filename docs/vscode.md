# Connecting to Biowulf using VScode

## Overview
Instructions adapted from the Biowulf tutorial [here](https://hpc.nih.gov/apps/vscode.html)
## Installing VScode

**Windows** users can install yourself from [here](https://code.visualstudio.com/docs/setup/windows).

**Mac** users should have IT install it for you from [here](https://code.visualstudio.com/docs/setup/mac).

## VScode setup
After installing `VScode`, the steps for setting it up to interface with Biowulf are

### Install the `Remote - SSH` extension in `VScode`
Click the `Extensions` icon in the left pane (four square blocks), search for the package, and install it.

### Ensuring `ssh-agent` is running
Windows: have IT start an administrative `PowerShell` and execute the following to allow you to start `ssh-agent`:
```
Get-Service -Name ssh-agent | Set-Service -StartupType Manual
```
Close the administrative `PowerShell` window when finished.

### Generate `ssh` private and public key pair on your local machine (laptop)
Start up `terminal` (macOS) or `PowerShell` (Windows) 
```bash
ssh-agent                   # start ssh-agent
ssh-keygen -t rsa -b 4096   # generate key pair
ssh-add                     # add the private key to your ssh-agent
```
You can press enter to continue through any prompts to automatically generate the keys `id_rsa` and `id_rsa.pub`. The public keys is the one ending with `.pub`.

### Adding the public key to your Biowulf `authorized_keys` file
First, copy your public key to a temporary file on Biowulf

**Windows PowerShell**: (be sure to replace `USERNAME` with your laptop username, and `COMPUTEID` with your biowulf compute ID):
```
scp C:/Users/USERNAME/.ssh/id_rsa.pub COMPUTEID@biowulf.nih.gov:~/key.tmp
```

**Mac Terminal**: (be sure to replace `COMPUTEID` with your Biowulf compute ID):
```
scp  ~/.ssh/id_rsa.pub COMPUTEID@biowulf.nih.gov:~/key.tmp
```

### Setting up the `ssh config` file on your laptop (and telling VScode where it is)
**Windows**: save the following as `C:\Users\USERNAME\.ssh\config`
```
Host cn*
User COMPUTEID
ProxyCommand C:\Windows\System32\OpenSSH\ssh.exe -o ForwardAgent=yes COMPUTEID@biowulf.nih.gov nc -w 120ms %h %p
```

**Mac**: save the following as `~/.ssh/config`
```
Host cn*
User COMPUTEID
ProxyCommand /usr/bin/ssh -o ForwardAgent=yes COMPUTEID@biowulf.nih.gov nc -w 120ms %h %p
```

## Connecting to Biowulf
Once `VScode` is configured, you will be able to connect to a *currently running* interactive session.
- Start an interactive session on Biowulf using `terminal` (macOS) or `MobaXterm` (Windows).
- Once you are granted your session, take note of the compute node name, something like `cn4244`.
- Within `VScode`, click the blue `><` icon in the lower left, then click `Connect to Host...`
- Enter the compute node ID of your interactive session, e.g. `cn4244`
- Click yes/accept whenever prompted. You should connect to Biowulf without needing to enter your password.
- If you are prompted for your password, your `ssh` keys are not properly configured.

## Editing scripts
The simplest way is to mount your `/data` drive following instructions here https://hpc.nih.gov/docs/helixdrive.html and then editing files using VScode.

[VScode web link](https://vscode.dev/)


## Setting up SSH keys

## 

