clear all;

% elgamal暗号化および復号化のための関数を定義

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

% 公開鍵と秘密鍵を生成
function [pk, sk] = elg_keygen(bits)
    % 素数p,qの生成
    while true
        q = nextprime(randi([2^(bits-2), 2^(bits-1)-1]));
        p = 2 * q + 1;
        if isprime(p)
            break;
        end
    end

    % 原始根gの生成
    while true
        g = randi([3, p-1]);
        if extpower(g, 2, p) ~= 1 && extpower(g, q, p) ~= 1
            break;
        end
    end
    % 秘密値xの生成
    x = randi([2, p-2]);

     % 公開値yの計算
    y = extpower(g, x, p);

    % 公開鍵と秘密鍵を返す
    pk = [p, g, y];
    sk = x;
end

% RSA暗号化関数
function enc = elg_encrypt(m, pk)
    % 行列mの要素を整数に変換してから暗号化
    gamma = 10;
    m = fix(m*gamma);

    p = pk(1);
    g = pk(2);
    y = pk(3);
    % ランダムなrを生成
    r = randi([2, p-2]);
    % 暗号文を計算(鍵情報のみ)
    % c1とc2は互いの値を参照しあわない
    c1 = extpower(g, r, p);

    e = extpower(y, r, p);
    if(m < 0)
        m = mod(m + p,p);
        c2 = -mod(m*e,p);
        enc = {c1,c2};
        return
    end
    c2 = mod(m*e,p);
    enc = {c1,c2};
end

% RSA復号化関数
function m = elg_decrypt(enc, pk, sk)
    gamma = 10;
    p = pk(1);
    c1 = enc{1};
    c2 = enc{2};
    
    if(c2 < 0)
        c2 = mod(c2 + p,p);
        e = extpower(c1, p - 1 - sk, p);
        m = -mod(c2*e,p)/gamma;
        return
    end
    e = extpower(c1, p - 1 - sk, p);
    m = mod(c2*e,p)/gamma;
end

function encenc = elg_times(enc1,enc2,pk,sk)
    gamma = 10;
    p = pk(1);

    encenc{1} = mod(enc1{1}*enc2{1},p);
    encenc{2} = mod(enc1{2}*enc2{2},p);

   encenc = elg_encrypt(elg_decrypt(encenc, pk, sk)/gamma,pk);
end

% 初期値応答を計算する関数を定義
function dxdt = state_eq(t, x, A, b, f)
    % 状態方程式 dx/dt = (A - b*f) * x
    dxdt = (A - b * f) * x;
end
% ===== メインプログラム =====

clear all; close all;

% ELGAMALキー生成
[pk, sk] = elg_keygen(18);
t = 0:0.01:4;

% システムパラメータを与える
A = [2 0.03; -0.01 0.1]; % 行列A
b = [1; 1]; % ベクトルb
c = [1 0; 0 1]; % x_1とx_2をプロットするためにcを単位行列にする

% システム行列Aの固有値を求める
eig(A); % システム行列Aの固有値

% 配置する極を与える
p1 = [-3; -3]; % -2, -3の場合
p2 = [-5; -6]; % -5, -6の場合

% アッカーマンの方法によりフィードバックゲインを求める
f1 = acker(A, b, p1); % -2, -3の場合のフィードバックゲインf1
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

% A_enc を 2x2 のセル配列として定義
A_enc = cell(2, 2);
% A行列の各要素を暗号化
for i = 1:2
    for j = 1:2
        % elg_encrypt が返す 1x2 のセル配列をそのまま代入
        A_enc{i, j} = elg_encrypt(A(i, j), pk);
    end
end

% A_dec を 2x2 のセル配列として定義
A_dec = zeros(2, 2);
for i =1:2
    for j = 1:2
        enc = A_enc(i,j);
        enc = enc{1};
        A_dec(i,j) = elg_decrypt(enc, pk,sk);
    end
end

% 状態変数の時系列データを暗号化
x1_enc = cell(length(t), 2);
x2_enc = cell(length(t), 2);
for j = 1:2
    for i = 1:length(t)
        x1_enc{i, j} = elg_encrypt(x1(i, j), pk);
        x2_enc{i, j} = elg_encrypt(x2(i, j), pk);
    end
end

% 状態変数の時系列データを復号化
x1_dec = zeros(length(t), 2);
x2_dec = zeros(length(t), 2);
for j = 1:2
    for i = 1:length(t)
        x1_dec(i, j) = elg_decrypt(x1_enc{i, j}, pk, sk);
        x2_dec(i, j) = elg_decrypt(x2_enc{i, j}, pk, sk);
    end
end

% フィードバックゲインの暗号化と復号化
f1_enc = cell(1, 2);
f2_enc = cell(1, 2);
f1_dec = zeros(1, 2);
f2_dec = zeros(1, 2);
for i = 1:2
    f1_enc{i} = elg_encrypt(f1(i), pk);
    f2_enc{i} = elg_encrypt(f2(i), pk);
end
for i = 1:2
    f1_dec(i) = elg_decrypt(f1_enc{i}, pk, sk);
    f2_dec(i) = elg_decrypt(f2_enc{i}, pk, sk);
end


for i = 1:length(t)
    % f1_enc と x1_enc の暗号化された値を掛け合わせ、復号化する
    f1x1_1_enc = elg_times(f1_enc{1}, x1_enc{i, 1}, pk, sk);
    f1x1_2_enc = elg_times(f1_enc{2}, x1_enc{i, 2}, pk, sk);

    % f2_enc と x2_enc の暗号化された値を掛け合わせ、復号化する
    f2x2_1_enc = elg_times(f2_enc{1}, x2_enc{i, 1}, pk, sk);
    f2x2_2_enc = elg_times(f2_enc{2}, x2_enc{i, 2}, pk, sk);
    
    % f1_enc と x1_enc の暗号化された値を足し合わせた後に復号化
    u1_dec(i) = -(elg_decrypt(f1x1_1_enc, pk, sk) + elg_decrypt(f1x1_2_enc, pk, sk));

    % f2_enc と x2_enc の暗号化された値を足し合わせた後に復号化
    u2_dec(i) = -(elg_decrypt(f2x2_1_enc, pk, sk) + elg_decrypt(f2x2_2_enc, pk, sk));
end

% y_decの初期化
y_dec = zeros(2, length(t));
y_dec(:,1) = x0;

% 状態方程式の計算（復号化された状態変数を使用して）
for i = 1:length(t)-1
    y_dec(:,i+1) = y_dec(:,i) + (A*x1_dec(i,:)' + b*u1_dec(i))*0.01;
end

y_enc = cell(1, length(t)); % セル配列に変更
plot_y_enc = zeros(1, length(t)); % プロット用の配列
% 暗号化とプロット用データの準備
for i = 1:length(t)
    y_enc{i} = elg_encrypt(y_dec(1,i), pk);
    plot_y_enc(i) = y_enc{i}{2}; % 各セルの第二要素を抽出
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
plot(t, plot_y_enc , '-b');
grid on; % 罫線を表示
xlabel('time t[s]'); % 横軸のラベル表示
ylabel('Enc_y'); % 縦軸のラベル表示
legend('y'); % 凡例の表示