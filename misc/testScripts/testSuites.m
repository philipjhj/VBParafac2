import matlab.unittest.TestSuite

%results = runtests(strcat('/../..',pwd),'Recursively',true);
%results

testDir = 'tests/';

myTests = TestSuite.fromFile(strcat(testDir,'qDistributionTest.m'));


%%
run(myTests)