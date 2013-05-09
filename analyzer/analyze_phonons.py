#!/usr/bin/env python

"""
This python script is a simple user interface to 
the C++ analyzer program "analyze_phonons.cpp".
To execute this script, the above c++ program has to be
compiled and made executable.
"""
import sys, os
import optparse
import datetime
import subprocess

parser = optparse.OptionParser()
parser.add_option('-t', '--temp')
parser.add_option('-m', '--mode')
parser.add_option('-k', '--kpoint')
parser.add_option('-c', '--calc')
parser.add_option('-l', '--length')
parser.add_option('-d', '--direction')

options, args = parser.parse_args()
file_result = args[0]

dir_obj = os.path.dirname(__file__)
analyze_obj = dir_obj + "/analyze_phonons "

if __name__ == '__main__':
    
    print "#", datetime.datetime.today()
    print "# Input file : ", file_result
    print "#"

    calc = options.calc

    if calc == "tau":
        if options.temp == None: # Damping function 
            if options.kpoint == None or options.mode == None:
                sys.exit("Invalid combination of --temp, --kpoint, and --mode options")
            else:
                if len(options.kpoint.split(':')) != 1:
                    sys.exit("Invalid usage of --kpoint for --calc=tau")
                if len(options.mode.split(':')) != 1:
                    sys.exit("Invalid usage of --mode for --calc=tau")

                target_k = int(options.kpoint)
                target_s = int(options.mode)
                calc = "tau_temp"
                command = analyze_obj + file_result + " " + calc + " " + str(target_k) + " " + str(target_s) 

                subprocess.call(command, shell=True)

        else: # Relaxation time (omega,tau)

            if options.kpoint == None:
                beg_k = 1
                end_k = 0
            else:
                if len(options.kpoint.split(':')) == 1:
                    beg_k = int(options.kpoint) 
                    end_k = beg_k 
                elif len(options.kpoint.split(':')) == 2:
                    arr = options.kpoint.split(':')
                    beg_k, end_k = int(arr[0]), int(arr[1])
                else:
                    sys.exit("Invalid usage of --kpoint for --calc=tau")

            if options.mode == None:
                beg_s = 1
                end_s = 0
            else:
                if len(options.mode.split(':')) == 1:
                    beg_s = int(options.mode)
                    end_s = beg_s
                elif len(options.mode.split(':')) == 2:
                    arr = options.mode.split(':')
                    beg_s, end_s = int(arr[0]), int(arr[1])
                else:
                    sys.exit("Invalid usage of --mode for --calc=tau")

            command = analyze_obj + file_result + " " + calc + " " + str(beg_k) + " " + str(end_k) \
                      + " " + str(beg_s) + " " + str(end_s) + " " + options.temp

            subprocess.call(command, shell=True)

    elif calc == "kappa_size":
        if options.temp == None:
            sys.exit("--temp should be given when --value = kappa_size")

        
        if not (options.kpoint == None):
            print "# Warning: --kpoint option is discarded"
        
        if options.mode == None:
            beg_s = 1
            end_s = 0
        else:
            if len(options.mode.split(':')) == 1:
                beg_s = int(options.mode) 
                end_s = beg_s
            elif len(options.mode.split(':')) == 2:
                arr = options.mode.split(':')
                beg_s, end_s = int(arr[0]), int(arr[1]) 
            else:
                sys.exit("Invalid usage of --mode for --calc=kappa_size")

        if options.length == None:
            max_len = 1000.0
            d_len = 1.0
        elif len(options.length.split(':')) == 2:
            arr = options.length.split(':')
            max_len, d_len = float(arr[0]), float(arr[1])
        else:
            sys.exit("Invalid usage of --length option")
            
        size_flag = [0]*3

        if options.direction == None:
            for i in range(3):
                size_flag[i] = 0
        else:
            if len(options.direction.split(':')) > 3:
                sys.exit("Invalid usage of --direction")

            arr = options.direction.split(':')
            for i in range(len(options.direction.split(':'))):
                size_flag[int(arr[i])-1] = 1

        command = analyze_obj + file_result + " " + calc + " " + str(beg_s) + " " + str(end_s) \
                  + " " + str(max_len) + " " + str(d_len) + " " + options.temp \
                  + " " + str(size_flag[0]) + " " + str(size_flag[1]) + " " + str(size_flag[2])
        subprocess.call(command, shell=True)

    elif calc == "kappa":

        if not (options.kpoint == None):
            print "# Warning: --kpoint option is discarded"

        if options.mode == None:
            beg_s = 1
            end_s = 0
        else:
            if len(options.mode.split(':')) == 1:
                beg_s = int(options.mode)
                end_s = beg_s 
            elif len(options.mode.split(':')) == 2:
                arr = options.mode.split(':')
                beg_s, end_s = int(arr[0]), int(arr[1]) 
            else:
                sys.exit("Invalid usage of --mode for --calc=kappa")

        command = analyze_obj + file_result + " " + calc + " " + str(beg_s) + " " + str(end_s)

        subprocess.call(command, shell=True)

    else:
        sys.exit("Invalid --calc option given")
