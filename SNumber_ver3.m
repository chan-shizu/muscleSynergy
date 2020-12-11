clc;
clear;
% pathOfMuscleActivity = "measuredWave\SwalOptResult_27m_ver3.csv";%�C���|�[�g����csv�t�@�C���̒���10^-18���炢��菬���Ȑ��́Amatlab�ł͌�����l�ɂȂ�̂Ŏ��O��csv�t�@�C�����C��
pathOfMuscleActivity = "measuredWave\arranged\arranged_Active_level�O�Սl�@63m.csv";
Wav = readmatrix(pathOfMuscleActivity);
DataNumW = size(Wav);
plot(Wav(:,:));

%% ���������̐ݒ�
d = 25; % �ؓd�Z���T�[�̐�
L = 15;%15; % �J��Ԃ��v�Z�̉�
dawn = 1;
dt = 0.01*dawn;%0.001*dawn; % ���Ԃ̍�
maxSynergyNumber = 10;
% movetime=3.0;%2.2;
movetime=DataNumW(1)*dt;%2.2;
%     caltime=3.5;
caltime=movetime;%2.2;
M2=0;
calculatedWave=0;
% Tw = 1.2/(dawn*dt);%0.9*1000/dawn ;   % �V�i�W�[�p������Tw�̐ݒ�
TwTemp = int8(caltime/(dawn*dt)/2);%0.9*1000/dawn ;   % �V�i�W�[�p������Tw�̐ݒ�
% Tw = double(Tw);
Tw = zeros(1,maxSynergyNumber);
Tw(:) = TwTemp;

startt=1; %����Ɍv�����J�n����Ă���_?
tmaxW = DataNumW(1);
FW=0;
FG=0;
endlw = tmaxW;

%% �֐��ŏ����g�`ST�쐬�@�����̕�����
makeST_ver2(Tw,dawn, maxSynergyNumber, pathOfMuscleActivity,d);

%% �f�[�^�̊Ԉ����C������
WavData = Wav;%(:,2:2:d*2);
WavData = movmean(WavData,10,1); %�ړ�����
WavData = downsample(WavData,dawn); %�f�[�^���Ԉ���
plot(WavData(:,:));

WavData = movmean(WavData,10,1); %�ړ�����
DataNumW = size(WavData);
tmaxW = DataNumW(1);
%     MA=max(max(WavData(1:movetime:1000/dawn,1:d)));%WavData�̒��ň�ԑ傫�Ȓl��MA�ɑ��
MA=max(max(WavData(:,:))); %�ÒJ�ύX�@WavData�̒��ň�ԑ傫�����Ŋ���ׂ��ł�?
WavData(:,:)=WavData(:,:)/MA;%���K��
WavDataT = WavData.';%�]�u
E2min = 10000; %?
rmax = 0 ;

for N=1:maxSynergyNumber
    N
    topt = zeros(1,2,N,L);
    calculatedWave=0;
    E2min = 10000; %?
    %     dt = 0.001*dawn; % ���ԍ�
    
    %% ���łȂ������g�`wi�̐ݒ�
    
    ST_name=['ST\ST',num2str(N)];
    ST_name_i =[ST_name,'.csv'];
    G1 = csvread(ST_name_i,0,0);
    
    w = zeros(d,max(Tw(1:N)),N);
    for i=1:N
        w(:,1:Tw(i),i) = G1(((i-1)*d+1):((i-1)*d+d),1:Tw(i)); %Tw(i) = 2.0*1000/dawn; % �V�i�W�[�p������Tw�̐ݒ�
    end
    
    w(w==0) = 10e-100;
    M=0;
    WH=0;
    
    c(:)=mean(mean(WavData))/mean(mean(w(:,:,1:N)));  %�W��ci�̌v�Z, �e�ؓ���(�ؓd�}�g�`�̕���/�v�Z�g�`�̕���)�̕���
%         c = ones(1,N);
    %% 1��ڂ̌v�Z
    for  l=1:L    %�J��Ԃ��v�Z�̉�
        c(i)=mean(mean(WavData))/mean(mean(w(:,:,i)));
        if l==1
            for i=1:N   % wnew2��3�����̐�
                
                zc=0;
                %             c(i)=mean(mean(WavData))/mean(mean(w(:,:,i))) ;  %�W��ci�̌v�Z, �e�ؓ���(�ؓd�}�g�`�̕���/�v�Z�g�`�̕���)�̕���
                
                if i==1              %���ԃV�t�gti�͈͎̔w��
                    synergyStartTime(i)=1; %�ŏ��̃V�i�W�[�͋؊����x�Ɠ����^�C�~���O�ŃX�^�[�g
                elseif i==N
                    synergyStartTime(i)=tmaxW-Tw(i);
                else
                    %                     if topt(1,1,i,l-1)+20+Tw(i)>tmaxW  %topt(1,1,i-1,l)+20+Tw(i)>tmaxW
                    %                         %                     synergyStartTime(i)=tmaxW-Tw(i);
                    %                         synergyStartTime(i)=tmaxW-Tw(i)-20;
                    %                     else
                    %                         %                     synergyStartTime(i)=topt(1,1,i-1,l)+20; %���̌�ŋ��߂�topt��20�𑫂�
                    synergyStartTime(i)=20; %���̌�ŋ��߂�topt��20�𑫂�
                    %               synergyStartTime(i)=topt(1,1,i-1)+5; %���̌�ŋ��߂�topt��20�𑫂�
                    %                     end
                end
                
                if i==1
                    %                 synergyEndTime(i)=0.1/dt;
                    synergyEndTime(i)=1;
                    
                elseif i==N
                    %                 if topt(1,1,i-1)+1*Tw(i)>tmaxW %�ؓd�}�̌v�����Ԃ�topt(1,1,i-1)+2*Tw��������Ƃ�
                    synergyEndTime(i)=tmaxW-Tw(i); %�Ō�̃V�i�W�[�̏I������
                    %                 else
                    %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    %                 end
                    
                else %1<i<N�̏ꍇ
                    
                    %                     if topt(1,1,i-1,l)+Tw(i)>tmaxW
                    %                         synergyEndTime(i)=tmaxW-Tw(i);
                    %                     else
                    %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    synergyEndTime(i)=Tw(i);
                    %                     end
                end
                
                for x=synergyStartTime(i):synergyEndTime(i)
                    
                    Matrix(1:d,1:tmaxW,i) = 0;
                    Matrix(1:d,x:Tw(i)+x-1,i) = w(:,1:Tw(i),i); %�V�i�W�[�p�����Ԃ�w�̒l��matrix�ɑ��
                    Plus=0;
                    
                    for y=1:tmaxW
                        Times = WavData(y,1:d)*Matrix(1:d,y,i); %�ϑ��g�`�Ɖ��̌v�Z�g�`�̊|���Z. �Ȃ�1:4? ���ׂĂ̋ؓ��Ȃ�1:7�ł���ׂ�. ���ԑ���
                        Plus = Plus+Times;
                    end
                    
                    Answer(x, (i-1)*d+1) = Plus;
                    
                    if Answer(x,(i-1)*d+1) > topt(1,2,i,l)
                        topt(1,1:2,i,l) = [x-1,Answer(x,(i-1)*d+1)];
                    end
                end
                
                Answer(1,(i-1)*d+2:(i-1)*d+3) = topt(1,1:2,i);         %���ԃV�t�gti�̌���
                WH(1:d,1:tmaxW,i) = 0;
                WH(1:d,topt(1,1,i)+1:topt(1,1,i)+Tw(i),i) =  c(i) * w(:,1:Tw(i),i); %���̌v���g�`
                M =   M+WH(:,:,i); %M�̍X�V
            end
        else
            %% topt�̌���
            for i=1:N   % wnew2��3�����̐�
                
                zc=0;
                %             c(i)=mean(mean(WavData))/mean(mean(w(:,:,i))) ;  %�W��ci�̌v�Z, �e�ؓ���(�ؓd�}�g�`�̕���/�v�Z�g�`�̕���)�̕���
                
                if i==1              %���ԃV�t�gti�͈͎̔w��
                    synergyStartTime(i)=1; %�ŏ��̃V�i�W�[�͋؊����x�Ɠ����^�C�~���O�ŃX�^�[�g
                elseif i==N
                    synergyStartTime(i)=tmaxW-Tw(i);
                else
                    %                     if topt(1,1,i-1,l)+20+Tw(i)>tmaxW  %topt(1,1,i-1,l)+20+Tw(i)>tmaxW
                    %                         %                     synergyStartTime(i)=tmaxW-Tw(i);
                    %                         synergyStartTime(i)=tmaxW-Tw(i);
                    if topt(1,1,i-1,l)>Tw(i)  %topt(1,1,i-1,l)+stepWidth+Tw(i)>tmaxW
                        synergyStartTime(i)=topt(1,1,i+1,l)-Tw(i);
                    else
                        %                     synergyStartTime(i)=topt(1,1,i-1,l)+20; %���̌�ŋ��߂�topt��20�𑫂�
                        synergyStartTime(i)=topt(1,1,i-1,l)+20; %���̌�ŋ��߂�topt��20�𑫂�
                        %               synergyStartTime(i)=topt(1,1,i-1)+5; %���̌�ŋ��߂�topt��20�𑫂�
                    end
                end
                
                if i==1
                    %                 synergyEndTime(i)=0.1/dt;
                    synergyEndTime(i)=1;
                    
                elseif i==N
                    %                 if topt(1,1,i-1)+1*Tw(i)>tmaxW %�ؓd�}�̌v�����Ԃ�topt(1,1,i-1)+2*Tw��������Ƃ�
                    synergyEndTime(i)=tmaxW-Tw(i); %�Ō�̃V�i�W�[�̏I������
                    %                 else
                    %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    %                 end
                    
                else %1<i<N�̏ꍇ
                    
                    if topt(1,1,i-1,l)+2*Tw(i)>tmaxW
                        synergyEndTime(i)=tmaxW-Tw(i);

                    else
                        %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                        synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    end
                end
                
                for x=synergyStartTime(i):synergyEndTime(i)
                    
                    Matrix(1:d,1:tmaxW,i) = 0;
                    Matrix(1:d,x:Tw(i)+x-1,i) = w(:,1:Tw(i),i); %�V�i�W�[�p�����Ԃ�w�̒l��matrix�ɑ��
                    Plus=0;
                    
                    % N�̃V�i�W�[ ���^����ꂽ�Ƃ��������̉��ŁC�Ƃ肤��S�Ă̎��ԃV�t�gTsi�ł�i�Ԗڂ̃V�i�W�[ ��s���s�ڂ̋؊����g�`ms(t) �Ƃ̎��ԑ��փ�si ���v�Z����D
                    for y=1:tmaxW
                        Times = WavData(y,1:d)*Matrix(1:d,y,i); %�ϑ��g�`�Ɖ��̌v�Z�g�`�̊|���Z. �Ȃ�1:4? ���ׂĂ̋ؓ��Ȃ�1:7�ł���ׂ�. ���ԑ���
                        Plus = Plus+Times;
                    end
                    
                    Answer(x, (i-1)*d+1) = Plus;
                    
                    %��A�̑����S�g���ׂĂ̊ϑ��g�`�ɑ΂��ČJ��Ԃ��C���ԃV�t�g������o���D
                    if Answer(x,(i-1)*d+1) > topt(1,2,i,l)
                        topt(1,1:2,i,l) = [x-1,Answer(x,(i-1)*d+1)];
                    end
                end
                
                Answer(1,(i-1)*d+2:(i-1)*d+3) = topt(1,1:2,i);         %���ԃV�t�gti�̌���
                WH(1:d,1:tmaxW,i) = 0;
                WH(1:d,topt(1,1,i)+1:topt(1,1,i)+Tw(i),i) =  c(i) * w(:,1:Tw(i),i); %���̌v���g�`
                M =   M+WH(:,:,i); %M�̍X�V
            end
        end
        
        %% �W��c�̌���
%                 for synergyNumber = 1:N
%         %                 updateC(synergyNumber)=(trace(WavData * M))./(trace(M.' * M)); %�f�R���|�W�V�����A���S���Y���ł̎c��
%                         updateC(synergyNumber) = mean(mean(WavData))/mean(mean(M));
%         %                 c(1:N) = c(1:N).*updateC ;
%                         c(synergyNumber) = c(synergyNumber)*updateC(synergyNumber) ;
%         %                   c(1:N) = cStart;
%         %                 c = cnew;
%                 end
        
        %% w�̍X�V
        for i=1:N
            wt(1:d,1:tmaxW,i) = 0;
            wt(1:d,topt(1,1,i)+1:topt(1,1,i)+Tw(i),i) = w(:,1:Tw(i),i);
            for x = 1:tmaxW
                if WH(1,x,i) == 0
                    wnew(1:d,x,i) = 0;
                else
                    wnew(1:d,x,i) = wt(:,x,i).*(WavDataT(:,x)./M(:,x)); %�؃V�i�W�[���X�V
                end
            end
            w(:,1:Tw(i),i) = wnew(:,topt(1,1,i,l)+1:topt(1,1,i,l)+Tw(i),i);
        end
        
        E2(l) = trace((WavDataT-calculatedWave).' * (WavDataT-calculatedWave)); %E^2 �f�R���|�W�V�����A���S���Y��
        
        % �J��Ԃ��v�Z�̂Ȃ���E2���ŏ��ɂȂ�l���L�^
        if (l>10) & (E2min > E2(l))
            E2min = E2(l);
            minNumber = l;
%             l
        end
        
        for i = 1:N
            wnew2Result((l-1)*d*N+(i-1)*d+1:(l-1)*d*N+(i-1)*d+d,:) = w(:,:,i) ;
            cResult((l-1)*(N+1)+i) = c(i); %�J��Ԃ��v�Z�̂��тɃV�i�W�[���̐��������������Ă����؊����x�̏I������
        end
    end
    %     Figure1 = figure();
    %     plot(wnew(5,:,i));
    %     Figure2 = figure();
    %     plot(M(5,:));
    
    csvwrite('Result.csv', Answer);
    csvwrite('w.csv', wnew2Result);
    csvwrite('c.csv', cResult);
    
    minNumber = 1;
    cend(1:N) = cResult((minNumber-1)*(N+1)+1:(minNumber-1)*(N+1)+N) ;%�e�V�i�W�[�̏I������, min�͌J��Ԃ��v�Z�̒��ōł����x������������
    wend(1:d*N,:) = wnew2Result((minNumber-1)*d*N+1:((minNumber-1)*d*N+d*N),:); %�ł�E2�������������Ƃ���wnew2Result�̒l��wend�Ƃ���
    tend(1:N) = topt(1,1,1:N,minNumber);
    
    %     synergyWave(1:(d+1)*N,1:tmaxW) = 0;
    synergyWave(1:d*N,1:tmaxW) = 0;
    
    for n=1:N % Tw�̕����������Ȃ�wend��tmaxW�̕��ŕ\��
        synergyWave(((n-1)*d+1):((n-1)*d+d),tend(1,n)+1:tend(1,n)+Tw(i)) = cend(1,n)* wend(((n-1)*d+1):((n-1)*d+d),1:Tw(i));
    end
    
    mend(d,tmaxW) = 0;
    
    for n=1:N %���ꂼ��̃V�i�W�[�̒l���ؓ����Ƃɑ������킹��
        mend = mend + synergyWave((1+d*(n-1)):(d*(n-1)+d),1:tmaxW);
    end
    
    meanMend = mean(mend(:,:),2);
    % ����W���̌v�Z
    %     V(:,:)=(1/(tmaxW))*(WavDataT(:,:)-Meanw2(:,:)).*(WavDataT(:,:)-Meanw2(:,:));
    %     V2 = sum(V,2); %�e�ؓ��̋ؓd�}�g�`�̕��U
    %     V3 = mean(V2); %�e�ؓ��̋ؓd�}�g�`�̕��U�̕���
    %
    %     W(:,:)=(1/(tmaxW))*(mend(:,:)-Mean2(:,:)).*(mend(:,:)-Mean2(:,:));
    %     W2 = sum(W,2); %�e�ؓ��̌v�Z�g�`�̕��U
    %     W3 = mean(W2); %�e�ؓ��̌v�Z�g�`�̕��U�̕���
    %     S(:,:) = W3 / V3 ;
    
    % ����W��(�l�̋N���ɂ����铮��v���~�e�B�u�̖���)
    sumOfDifference = mean((WavData(:,:) - mend(:,:).').^2); %�v���g�`�ƌv�Z�g�`�̌덷�̓��a
    V(:,:)=(1/(tmaxW))*(WavDataT(:,:)-meanMend(:,:)).*(WavDataT(:,:)-meanMend(:,:));
    V2 = sum(V,2); %�e�ؓ��̋ؓd�}�g�`�̕��U
    Seach(1:d) = 1-(sumOfDifference.^2)/((V2.').^2); %����W��
    S = min(Seach(1,:)); %�ł��덷���傫���ؓ���Seach�̒l
    
    csvwrite('synergyWave.csv', synergyWave); %�e�V�i�W�[�̊����x�𑫂����v�Z�g�`(�ؓ��̐�x�v�Z�X�e�b�v)
    SS(1,N)=N;
    SS(2,N)= S;
    
    % �V�i�W�[����1�������ꍇ�Ɣ�ׂĐ��x�����サ�Ă��Ȃ�������I��
    %     if (N>2) & (abs(SS(2,N)-SS(2,N-1)) < 0.05)
    %         if (max(max(synergyWave((N-2)*d+1:(N-2)*d+d,:)))-max(max(synergyWave((N-1)*d+1:(N-1)*d+d,:))) > 0.1)
    %             Snumber=SS(1,N);
    %         else
    %             if N>2
    %                 Snumber=SS(1,N-1);
    %             else
    %                 Snumber=SS(1,N);
    %             end
    %         end
    %         Snumber
    %         break
    %     end
    
    if N==3
        Snumber=SS(1,N);
        break
    end
    
    clearvars -except SS trials tlate selpath_i file_name file_name_i Wav WavData WavDataT MA d L Tw tmaxW tmax dt dawn caltime:
end

% �ϑ��g�`�̕\��
WavDataT = WavData';
plotType = 3
% make_plot(plotType,d,WavData,synergyWave,mend,N);
make_plot();


