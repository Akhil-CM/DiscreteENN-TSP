#!/usr/bin/python3

import os
import sys
import time
import argparse
import subprocess

COMMAND = ""

def runCommand(cmd, args):
    process = subprocess.Popen(
        [cmd] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    command_in_loop = 0
    while True:
        # if command_in_loop > 1e6:
        #     return 1
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
            # command_in_loop += 1
            if "[Error]" in output:
                return 1

    exit_code = process.wait()
    return exit_code

def runUntilFailure(cmd, args):
    count = 0
    while True:
        print()
        exit_code = runCommand(cmd, args)
        time.sleep(1)
        if exit_code != 0:
            print(f"\nExit code {exit_code} returned")
            break
        count += 1
        os.system('clear')
        print(f"\rRun: {count} successful", end='', flush=True)
        time.sleep(5)

    print(f"Completed runs: {count}")

def main():
    # parser = argparse.ArgumentParser(description="Run a command repeatedly until it fails.")

    # parser.add_argument("command_str",
    #                     type=str,
    #                     help="The command that should be run repeatedly.")
    # parser.add_argument("command_args", nargs='*', help="Additional arguments to pass to the command.")

    # args = parser.parse_args()

    # command = args.command_str
    # command_args = args.command_args

    # if not command:
    #     print(f"[Error]: Not command argument provided.\nExiting")
    #     exit(1)

    if len(sys.argv) < 2:
        print(f"[Error]: Not command argument provided.\nExiting")
        sys.exit(1)

    command = sys.argv[1]

    command_args = sys.argv[2:]

    if command[0:2] == "./":
        command = os.path.abspath(command)

    print()
    print()
    runUntilFailure(command, command_args)

if __name__ == "__main__":
    main()
