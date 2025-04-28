% run_all_tests_in_test_browser.m
% This script creates a test suite from your MATLAB library and
% opens the Test Browser app with the suite loaded.
%
% For details, see:
% https://www.mathworks.com/help/matlab/ref/testbrowser-app.html
% https://www.mathworks.com/help/matlab/ref/runtests.html
% https://www.mathworks.com/help/matlab/matlab_prog/run-tests-using-test-browser.html

clearvars; clc;

restoredefaultpath;

% Create a test suite from the current directory and its subfolders.
suite = testsuite(pwd, 'IncludeSubfolders', true);

% Launch the Test Browser app with the test suite loaded.
run(suite);

% Note:
% Currently, there is no public API to automatically trigger test execution
% within the Test Browser app. After the app loads your tests, you will need
% to click the "Run All Tests" button manually.
