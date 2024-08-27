clear all;

% RSA暗号化および復号化のための関数を定義


% 公開鍵と秘密鍵を生成
function [n, e, d] = rsa_keygen()
    p = 811; % 1つ目の素数
    q = 2027; % 2つ目の素数
    n = p * q; % nを計算
    l = lcm(p-1,q-1);
    e = 2;
    while(gcd(l,e) ~= 1)
        e = e+1;
    end
    d = 2;
    while(mod(e * d,l) ~= 1)
        d = d + 1;
    end
end

% 繰り返し自乗法を使った法nのべき乗計算（aのk乗をnで割った余りを求める）
function value = extpower(a,k,n)
    a = mod(a,n);
    if( (a == 0) || (n == 0))
        value = 0;
        return
    end
    if(k == 0)
        value = mod(1,n);
        return
    end
    value = 1;
    for i = 0:k-1
        value = value*a;
        if(value >= n)
            value = mod(value,n);
        end
    end
    return
end



% RSA暗号化関数
function c = rsa_encrypt(m, e, n)
    % 行列mの要素を整数に変換してから暗号化
    gamma = 100;
    m = fix(m*gamma);
    if(m < 0)
        m = mod(m + n,n);
        c = -extpower(m,e,n);
        return
    end
    
    c = extpower(m,e,n);
end

% RSA復号化関数
function m = rsa_decrypt(c, d, n)
    gamma = 100;
    if(c < 0)
        c = mod(c + n,n);
        m = -extpower(c,d,n)/gamma;
        return
    end
    m = extpower(c,d,n);
    m = m/gamma;
end

function ee = rsa_times(e1,e2,e,d,n)
    gamma = 100;
    ee = rsa_encrypt(rsa_decrypt(e1*e2,d,n)/gamma,e,n);
end

% 初期値応答を計算する関数を定義
function dxdt = state_eq(t, x, A, b, f)
    % 状態方程式 dx/dt = (A - b*f) * x
    dxdt = (A - b * f) * x;
end
% ===== メインプログラム =====

clear all; close all;

% RSAキー生成
[n, e, d] = rsa_keygen();
t = 0:0.01:4;

% システムパラメータを与える
A = [2 3; -1 4]; % 行列A
b = [1; 1]; % ベクトルb
c = [1 0; 0 1]; % x_1とx_2をプロットするためにcを単位行列にする

% システム行列Aの固有値を求める
eig(A); % システム行列Aの固有値

% 配置する極を与える
p1 = [-3; -3]; % -2, -3の場合
p2 = [-5; -6]; % -5, -6の場合

% アッカーマンの方法によりフィードバックゲインを求める
f1 = acker(A, b, p1) % -2, -3の場合のフィードバックゲインf1
f2 = acker(A, b, p2); % -5, -6の場合のフィードバックゲインf2

%%ここから

% システムの初期値を与える
x0 = [1; 1]; % 初期ベクトル


% 初期値応答の計算
[t, x1] = ode45(@(t, x) state_eq(t, x, A, b, f1), t, x0); % f1の場合の初期値応答
[t, x2] = ode45(@(t, x) state_eq(t, x, A, b, f2), t, x0); % f2の場合の初期値応答

% 制御入力の計算
u1 = -f1 * x1'; % f1の場合の制御入力
u2 = -f2 * x2'; % f2の場合の制御入力

% A行列の各要素を暗号化
for i =1:2
    for j = 1:2
        A_enc(i,j) = rsa_encrypt(A(i,j), e, n);
    end
end

for i =1:2
    for j = 1:2
        A_dec(i,j) = rsa_decrypt(A_enc(i,j), d, n);
    end
end

for j = 1:2
    for i=1:length(t)
        % 状態変数の時系列データを暗号化
        x1_enc(i,j)  = rsa_encrypt(x1(i,j) , e, n);
        x2_enc(i,j)  = rsa_encrypt(x2(i,j) , e, n);
    end
end
for j = 1:2
    for i=1:length(t)
        % 状態変数の時系列データを復号化
        x1_dec(i,j)  = rsa_decrypt(x1_enc(i,j) , d, n);
        x2_dec(i,j)  = rsa_decrypt(x2_enc(i,j) , d, n);
    end
end

for i = 1:2
    f1_enc(1,i) = rsa_encrypt(f1(1,i), e, n);
    f2_enc(1,i) = rsa_encrypt(f2(1,i), e, n);
end
for i = 1:2
    f1_dec(1,i) = rsa_decrypt(f1_enc(1,i), d, n);
    f2_dec(1,i) = rsa_decrypt(f2_enc(1,i), d, n);
end

% 制御入力の時系列データを暗号化
for i = 1:length(t)
    % 各要素enc*enc
    temp(i,1) = f1_enc(1,1)*x1_enc(i,1);
    temp(i,2) = f1_enc(1,2)*x1_enc(i,2);

    % 各要素enc*encしたあとにそのまま足し合わせる→これはできない．decした後じゃないと足せない．
    temp(i,3) = rsa_times(f1_enc(1,1),x1_enc(i,1),e,d,n)+rsa_times(f1_enc(1,2),x1_enc(i,2),e,d,n);
    temp(i,4) = rsa_times(f2_enc(1,1),x2_enc(i,1),e,d,n) + rsa_times(f2_enc(1,2),x2_enc(i,2),e,d,n);

    u1_enc = f1_enc(1,1)*x1_enc(i,1) + f1_enc(1,2)*x1_enc(i,2);
    u2_enc = f2_enc(1,1)*x2_enc(i,1) + f2_enc(1,2)*x2_enc(i,2);
end

for i=1:length(t)
    % 各要素enc*encしたあとにそれぞれdec
    temp(i,5) = rsa_decrypt(rsa_times(f1_enc(1,1),x1_enc(i,1),e,d,n),d,n); % 各要素enc*encをdec
    temp(i,6) = rsa_decrypt(rsa_times(f1_enc(1,2),x1_enc(i,2),e,d,n),d,n); % 各要素enc*encをdec

    %各要素のenc*encを足し合わせた後にdec→これはできない
    temp(i,9) = rsa_decrypt(rsa_times(f1_enc(1,1),x1_enc(i,1),e,d,n)+rsa_times(f1_enc(1,2),x1_enc(i,2),e,d,n),d,n);
    temp(i,10) = rsa_decrypt(rsa_times(f2_enc(1,1),x2_enc(i,1),e,d,n)+rsa_times(f2_enc(1,2),x2_enc(i,2),e,d,n),d,n);

    %各要素をenc*encしたあとににdecし，その後足し合わせ→これはできる
    u1_dec(1,i) = -(rsa_decrypt(rsa_times(f1_enc(1,1),x1_enc(i,1),e,d,n),d,n)+rsa_decrypt(rsa_times(f1_enc(1,2),x1_enc(i,2),e,d,n),d,n));
    u2_dec(1,i) = -(rsa_decrypt(rsa_times(f2_enc(1,1),x2_enc(i,1),e,d,n),d,n)+rsa_decrypt(rsa_times(f2_enc(1,2),x2_enc(i,2),e,d,n),d,n));
end

y_dec = zeros(2, length(t));
y_dec(:,1) = x0;
for i=1:length(t)-1
    y_dec(:,i+1) = y_dec(:,i) + (A*x1_dec(i,:)' + b*u1_dec(1,i))*0.01;
end

for i=1:length(t)
    y_enc(1,i) = rsa_encrypt(y_dec(1,i),e,n);
end

% 図9.1のプロット
figure(1) % 図のウィンドウを開く
plot(t, x1(:,1), '-b',t, y_dec(1,:), '-r'); % x_1をプロット
grid on; % 罫線を表示
xlabel('time t[s]'); % 横軸のラベル表示
ylabel('y(t)'); % 縦軸のラベル表示
legend('output without encryption','output y'); % 凡例の表示

% 図9.4のプロット（入力の大きさ）
figure(2) % 図のウィンドウを開く
plot(t, u1, '-b',t, u1_dec,'-r');
grid on; % 罫線を表示
xlabel('time t[s]'); % 横軸のラベル表示
ylabel('Control Input u(t)'); % 縦軸のラベル表示
legend('input without encryption','input u'); % 凡例の表示

% 図9.4のプロット（入力の大きさ）
figure(3) % 図のウィンドウを開く
plot(t, y_enc(1,:), '-b');
grid on; % 罫線を表示
xlabel('time t[s]'); % 横軸のラベル表示
ylabel('Enc_y'); % 縦軸のラベル表示
legend('y'); % 凡例の表示
