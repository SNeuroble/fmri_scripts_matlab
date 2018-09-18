% create clf imgs for all tasks
% before starting, run: sshfs -o IdentityFile=~/.ssh/MRCInstance1.pem
%    ec2-user@52.87.169.145:data/hcpTask/ mnt/

task={'SOCIAL'; 'WM'; 'GAMBLING'; 'RELATIONAL'; 'EMOTION'};
% task={'WM'};
for i=1:length(task)
    [tp_summary{i}]=counting_TPs(task{i});
%     [tps2{i}]=counting_TPs(task{i},tp_summary_new{i});
end