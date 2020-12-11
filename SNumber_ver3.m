clc;
clear;
% pathOfMuscleActivity = "measuredWave\SwalOptResult_27m_ver3.csv";%インポートしたcsvファイルの中の10^-18くらいより小さな数は、matlabでは誤った値になるので事前にcsvファイルを修正
pathOfMuscleActivity = "measuredWave\arranged\arranged_Active_level軌跡考察63m.csv";
Wav = readmatrix(pathOfMuscleActivity);
DataNumW = size(Wav);
plot(Wav(:,:));

%% 初期条件の設定
d = 25; % 筋電センサーの数
L = 15;%15; % 繰り返し計算の回数
dawn = 1;
dt = 0.01*dawn;%0.001*dawn; % 時間の差
maxSynergyNumber = 10;
% movetime=3.0;%2.2;
movetime=DataNumW(1)*dt;%2.2;
%     caltime=3.5;
caltime=movetime;%2.2;
M2=0;
calculatedWave=0;
% Tw = 1.2/(dawn*dt);%0.9*1000/dawn ;   % シナジー継続時間Twの設定
TwTemp = int8(caltime/(dawn*dt)/2);%0.9*1000/dawn ;   % シナジー継続時間Twの設定
% Tw = double(Tw);
Tw = zeros(1,maxSynergyNumber);
Tw(:) = TwTemp;

startt=1; %正常に計測が開始されている点?
tmaxW = DataNumW(1);
FW=0;
FG=0;
endlw = tmaxW;

%% 関数で初期波形ST作成　ただの放物線
makeST_ver2(Tw,dawn, maxSynergyNumber, pathOfMuscleActivity,d);

%% データの間引き，平滑化
WavData = Wav;%(:,2:2:d*2);
WavData = movmean(WavData,10,1); %移動平均
WavData = downsample(WavData,dawn); %データを間引き
plot(WavData(:,:));

WavData = movmean(WavData,10,1); %移動平均
DataNumW = size(WavData);
tmaxW = DataNumW(1);
%     MA=max(max(WavData(1:movetime:1000/dawn,1:d)));%WavDataの中で一番大きな値をMAに代入
MA=max(max(WavData(:,:))); %静谷変更　WavDataの中で一番大きい数で割るべきでは?
WavData(:,:)=WavData(:,:)/MA;%正規化
WavDataT = WavData.';%転置
E2min = 10000; %?
rmax = 0 ;

for N=1:maxSynergyNumber
    N
    topt = zeros(1,2,N,L);
    calculatedWave=0;
    E2min = 10000; %?
    %     dt = 0.001*dawn; % 時間差
    
    %% 負でない初期波形wiの設定
    
    ST_name=['ST\ST',num2str(N)];
    ST_name_i =[ST_name,'.csv'];
    G1 = csvread(ST_name_i,0,0);
    
    w = zeros(d,max(Tw(1:N)),N);
    for i=1:N
        w(:,1:Tw(i),i) = G1(((i-1)*d+1):((i-1)*d+d),1:Tw(i)); %Tw(i) = 2.0*1000/dawn; % シナジー継続時間Twの設定
    end
    
    w(w==0) = 10e-100;
    M=0;
    WH=0;
    
    c(:)=mean(mean(WavData))/mean(mean(w(:,:,1:N)));  %係数ciの計算, 各筋肉の(筋電図波形の平均/計算波形の平均)の平均
%         c = ones(1,N);
    %% 1回目の計算
    for  l=1:L    %繰り返し計算の回数
        c(i)=mean(mean(WavData))/mean(mean(w(:,:,i)));
        if l==1
            for i=1:N   % wnew2の3次元の数
                
                zc=0;
                %             c(i)=mean(mean(WavData))/mean(mean(w(:,:,i))) ;  %係数ciの計算, 各筋肉の(筋電図波形の平均/計算波形の平均)の平均
                
                if i==1              %時間シフトtiの範囲指定
                    synergyStartTime(i)=1; %最初のシナジーは筋活性度と同じタイミングでスタート
                elseif i==N
                    synergyStartTime(i)=tmaxW-Tw(i);
                else
                    %                     if topt(1,1,i,l-1)+20+Tw(i)>tmaxW  %topt(1,1,i-1,l)+20+Tw(i)>tmaxW
                    %                         %                     synergyStartTime(i)=tmaxW-Tw(i);
                    %                         synergyStartTime(i)=tmaxW-Tw(i)-20;
                    %                     else
                    %                         %                     synergyStartTime(i)=topt(1,1,i-1,l)+20; %この後で求めるtoptに20を足す
                    synergyStartTime(i)=20; %この後で求めるtoptに20を足す
                    %               synergyStartTime(i)=topt(1,1,i-1)+5; %この後で求めるtoptに20を足す
                    %                     end
                end
                
                if i==1
                    %                 synergyEndTime(i)=0.1/dt;
                    synergyEndTime(i)=1;
                    
                elseif i==N
                    %                 if topt(1,1,i-1)+1*Tw(i)>tmaxW %筋電図の計測時間をtopt(1,1,i-1)+2*Twが超えるとき
                    synergyEndTime(i)=tmaxW-Tw(i); %最後のシナジーの終了時間
                    %                 else
                    %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    %                 end
                    
                else %1<i<Nの場合
                    
                    %                     if topt(1,1,i-1,l)+Tw(i)>tmaxW
                    %                         synergyEndTime(i)=tmaxW-Tw(i);
                    %                     else
                    %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    synergyEndTime(i)=Tw(i);
                    %                     end
                end
                
                for x=synergyStartTime(i):synergyEndTime(i)
                    
                    Matrix(1:d,1:tmaxW,i) = 0;
                    Matrix(1:d,x:Tw(i)+x-1,i) = w(:,1:Tw(i),i); %シナジー継続時間のwの値をmatrixに代入
                    Plus=0;
                    
                    for y=1:tmaxW
                        Times = WavData(y,1:d)*Matrix(1:d,y,i); %観測波形と仮の計算波形の掛け算. なぜ1:4? すべての筋肉なら1:7であるべき. 時間相関
                        Plus = Plus+Times;
                    end
                    
                    Answer(x, (i-1)*d+1) = Plus;
                    
                    if Answer(x,(i-1)*d+1) > topt(1,2,i,l)
                        topt(1,1:2,i,l) = [x-1,Answer(x,(i-1)*d+1)];
                    end
                end
                
                Answer(1,(i-1)*d+2:(i-1)*d+3) = topt(1,1:2,i);         %時間シフトtiの決定
                WH(1:d,1:tmaxW,i) = 0;
                WH(1:d,topt(1,1,i)+1:topt(1,1,i)+Tw(i),i) =  c(i) * w(:,1:Tw(i),i); %仮の計測波形
                M =   M+WH(:,:,i); %Mの更新
            end
        else
            %% toptの決定
            for i=1:N   % wnew2の3次元の数
                
                zc=0;
                %             c(i)=mean(mean(WavData))/mean(mean(w(:,:,i))) ;  %係数ciの計算, 各筋肉の(筋電図波形の平均/計算波形の平均)の平均
                
                if i==1              %時間シフトtiの範囲指定
                    synergyStartTime(i)=1; %最初のシナジーは筋活性度と同じタイミングでスタート
                elseif i==N
                    synergyStartTime(i)=tmaxW-Tw(i);
                else
                    %                     if topt(1,1,i-1,l)+20+Tw(i)>tmaxW  %topt(1,1,i-1,l)+20+Tw(i)>tmaxW
                    %                         %                     synergyStartTime(i)=tmaxW-Tw(i);
                    %                         synergyStartTime(i)=tmaxW-Tw(i);
                    if topt(1,1,i-1,l)>Tw(i)  %topt(1,1,i-1,l)+stepWidth+Tw(i)>tmaxW
                        synergyStartTime(i)=topt(1,1,i+1,l)-Tw(i);
                    else
                        %                     synergyStartTime(i)=topt(1,1,i-1,l)+20; %この後で求めるtoptに20を足す
                        synergyStartTime(i)=topt(1,1,i-1,l)+20; %この後で求めるtoptに20を足す
                        %               synergyStartTime(i)=topt(1,1,i-1)+5; %この後で求めるtoptに20を足す
                    end
                end
                
                if i==1
                    %                 synergyEndTime(i)=0.1/dt;
                    synergyEndTime(i)=1;
                    
                elseif i==N
                    %                 if topt(1,1,i-1)+1*Tw(i)>tmaxW %筋電図の計測時間をtopt(1,1,i-1)+2*Twが超えるとき
                    synergyEndTime(i)=tmaxW-Tw(i); %最後のシナジーの終了時間
                    %                 else
                    %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    %                 end
                    
                else %1<i<Nの場合
                    
                    if topt(1,1,i-1,l)+2*Tw(i)>tmaxW
                        synergyEndTime(i)=tmaxW-Tw(i);

                    else
                        %                     synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                        synergyEndTime(i)=topt(1,1,i-1,l)+Tw(i);
                    end
                end
                
                for x=synergyStartTime(i):synergyEndTime(i)
                    
                    Matrix(1:d,1:tmaxW,i) = 0;
                    Matrix(1:d,x:Tw(i)+x-1,i) = w(:,1:Tw(i),i); %シナジー継続時間のwの値をmatrixに代入
                    Plus=0;
                    
                    % N個のシナジー が与えられたという条件の下で，とりうる全ての時間シフトTsiでのi番目のシナジー とs試行目の筋活動波形ms(t) との時間相関ψsi を計算する．
                    for y=1:tmaxW
                        Times = WavData(y,1:d)*Matrix(1:d,y,i); %観測波形と仮の計算波形の掛け算. なぜ1:4? すべての筋肉なら1:7であるべき. 時間相関
                        Plus = Plus+Times;
                    end
                    
                    Answer(x, (i-1)*d+1) = Plus;
                    
                    %一連の操作をS組すべての観測波形に対して繰り返し，時間シフトを割り出す．
                    if Answer(x,(i-1)*d+1) > topt(1,2,i,l)
                        topt(1,1:2,i,l) = [x-1,Answer(x,(i-1)*d+1)];
                    end
                end
                
                Answer(1,(i-1)*d+2:(i-1)*d+3) = topt(1,1:2,i);         %時間シフトtiの決定
                WH(1:d,1:tmaxW,i) = 0;
                WH(1:d,topt(1,1,i)+1:topt(1,1,i)+Tw(i),i) =  c(i) * w(:,1:Tw(i),i); %仮の計測波形
                M =   M+WH(:,:,i); %Mの更新
            end
        end
        
        %% 係数cの決定
%                 for synergyNumber = 1:N
%         %                 updateC(synergyNumber)=(trace(WavData * M))./(trace(M.' * M)); %デコンポジションアルゴリズムでの残差
%                         updateC(synergyNumber) = mean(mean(WavData))/mean(mean(M));
%         %                 c(1:N) = c(1:N).*updateC ;
%                         c(synergyNumber) = c(synergyNumber)*updateC(synergyNumber) ;
%         %                   c(1:N) = cStart;
%         %                 c = cnew;
%                 end
        
        %% wの更新
        for i=1:N
            wt(1:d,1:tmaxW,i) = 0;
            wt(1:d,topt(1,1,i)+1:topt(1,1,i)+Tw(i),i) = w(:,1:Tw(i),i);
            for x = 1:tmaxW
                if WH(1,x,i) == 0
                    wnew(1:d,x,i) = 0;
                else
                    wnew(1:d,x,i) = wt(:,x,i).*(WavDataT(:,x)./M(:,x)); %筋シナジーを更新
                end
            end
            w(:,1:Tw(i),i) = wnew(:,topt(1,1,i,l)+1:topt(1,1,i,l)+Tw(i),i);
        end
        
        E2(l) = trace((WavDataT-calculatedWave).' * (WavDataT-calculatedWave)); %E^2 デコンポジションアルゴリズム
        
        % 繰り返し計算のなかでE2が最小になる値を記録
        if (l>10) & (E2min > E2(l))
            E2min = E2(l);
            minNumber = l;
%             l
        end
        
        for i = 1:N
            wnew2Result((l-1)*d*N+(i-1)*d+1:(l-1)*d*N+(i-1)*d+d,:) = w(:,:,i) ;
            cResult((l-1)*(N+1)+i) = c(i); %繰り返し計算のたびにシナジー数の数だけ書き足していく筋活性度の終了時間
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
    cend(1:N) = cResult((minNumber-1)*(N+1)+1:(minNumber-1)*(N+1)+N) ;%各シナジーの終了時間, minは繰り返し計算の中で最も精度が高かった時
    wend(1:d*N,:) = wnew2Result((minNumber-1)*d*N+1:((minNumber-1)*d*N+d*N),:); %最もE2が小さかったときのwnew2Resultの値をwendとする
    tend(1:N) = topt(1,1,1:N,minNumber);
    
    %     synergyWave(1:(d+1)*N,1:tmaxW) = 0;
    synergyWave(1:d*N,1:tmaxW) = 0;
    
    for n=1:N % Twの幅しか持たないwendをtmaxWの幅で表現
        synergyWave(((n-1)*d+1):((n-1)*d+d),tend(1,n)+1:tend(1,n)+Tw(i)) = cend(1,n)* wend(((n-1)*d+1):((n-1)*d+d),1:Tw(i));
    end
    
    mend(d,tmaxW) = 0;
    
    for n=1:N %それぞれのシナジーの値を筋肉ごとに足し合わせる
        mend = mend + synergyWave((1+d*(n-1)):(d*(n-1)+d),1:tmaxW);
    end
    
    meanMend = mean(mend(:,:),2);
    % 決定係数の計算
    %     V(:,:)=(1/(tmaxW))*(WavDataT(:,:)-Meanw2(:,:)).*(WavDataT(:,:)-Meanw2(:,:));
    %     V2 = sum(V,2); %各筋肉の筋電図波形の分散
    %     V3 = mean(V2); %各筋肉の筋電図波形の分散の平均
    %
    %     W(:,:)=(1/(tmaxW))*(mend(:,:)-Mean2(:,:)).*(mend(:,:)-Mean2(:,:));
    %     W2 = sum(W,2); %各筋肉の計算波形の分散
    %     W3 = mean(W2); %各筋肉の計算波形の分散の平均
    %     S(:,:) = W3 / V3 ;
    
    % 決定係数(人の起立における動作プリミティブの役割)
    sumOfDifference = mean((WavData(:,:) - mend(:,:).').^2); %計測波形と計算波形の誤差の二乗和
    V(:,:)=(1/(tmaxW))*(WavDataT(:,:)-meanMend(:,:)).*(WavDataT(:,:)-meanMend(:,:));
    V2 = sum(V,2); %各筋肉の筋電図波形の分散
    Seach(1:d) = 1-(sumOfDifference.^2)/((V2.').^2); %決定係数
    S = min(Seach(1,:)); %最も誤差が大きい筋肉のSeachの値
    
    csvwrite('synergyWave.csv', synergyWave); %各シナジーの活性度を足した計算波形(筋肉の数x計算ステップ)
    SS(1,N)=N;
    SS(2,N)= S;
    
    % シナジー数が1小さい場合と比べて精度が向上していなかったら終了
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

% 観測波形の表示
WavDataT = WavData';
plotType = 3
% make_plot(plotType,d,WavData,synergyWave,mend,N);
make_plot();


