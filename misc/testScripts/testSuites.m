import matlab.unittest.TestSuite

%results = runtests(strcat('/../..',pwd),'Recursively',true);
%results    

testDir = 'tests/';

% myTests = TestSuite.fromFile(strcat(testDir,'qDistributionTest.m'));
myTests = TestSuite.fromClass(?qDistributionTest,'Tag','CAVI');
% myTests = TestSuite.fromClass(?qDistributionTest,'Tag','Updates');


results = run(myTests);

folderpath = 'output/tests/';
testname = 'all_updates';
save(strcat(folderpath,datestr(now,'yyyy_mm_dd_HHMM'),'_',testname));

load gong.mat;
for i = 1:3
    sound(y);
    pause(6)
end
