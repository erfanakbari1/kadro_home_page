% نمونه کد MATLAB برای محاسبه ضریب فوگاسیته، فاکتور تراکم‌پذیری و حجم مولی جزئی متان در مخلوط متان/اتان با استفاده از معادله‌ی حالت پاتل–تژا
% این کد تنها جهت راهنما است و ممکن است نیاز به اصلاح مقادیر ثابت‌ها داشته باشد.

clear; clc;

%% خواص بحرانی و عامل بی‌بعدی (acentric factor)
% [Tc (K), Pc (bar), omega]
comp(1).name  = 'Methane';
comp(1).Tc    = 190.6;        % K
comp(1).Pc    = 45.99;        % bar
comp(1).omega = 0.011;

comp(2).name  = 'Ethane';
comp(2).Tc    = 305.3;        % K
comp(2).Pc    = 48.84;        % bar
comp(2).omega = 0.099;

R = 0.08314;  % L·bar/(mol·K)

%% درصد مولی اجزا
z = [0.3, 0.7];

%% بازه‌های فشار (bar) و دما (K)
Tlist = linspace(280, 380, 5);   % مثال
Plist = linspace(10, 70, 5);

%% ثابت‌های معادله پاتل–تژا برای هر جزء
for i = 1:2
    Tc = comp(i).Tc;
    Pc = comp(i).Pc;
    omega = comp(i).omega;
    % پارامترهای اصلی
    comp(i).b = 0.07780 * R * Tc / Pc;
    kappa     = 0.37464 + 1.54226*omega - 0.26992*omega^2; % ممکن است نیاز به تغییر برای PT
    comp(i).a = 0.45724 * (R*Tc)^2 / Pc;
    comp(i).eps = 0.0 + 1.0*omega;   % تقریب؛ به مقالات PT مراجعه شود
    comp(i).sig = 1.0 - 1.0*omega;   % تقریب
    comp(i).kappa = kappa;
end

%% محاسبه برای هر دما و فشار
results = [];
for T = Tlist
    for P = Plist
        % پارامترهای a و b مخلوط
        a_i = zeros(1,2); b_i = zeros(1,2);
        eps_i = zeros(1,2); sig_i = zeros(1,2);
        for i=1:2
            alpha = (1 + comp(i).kappa*(1 - sqrt(T/comp(i).Tc)))^2;
            a_i(i) = comp(i).a * alpha;
            b_i(i) = comp(i).b;
            eps_i(i) = comp(i).eps;
            sig_i(i) = comp(i).sig;
        end
        % قوانین اختلاط ساده
        b_mix = sum(z .* b_i);
        a_mix = 0;
        for i=1:2
            for j=1:2
                a_mix = a_mix + z(i)*z(j)*sqrt(a_i(i)*a_i(j));
            end
        end
        eps_mix = sum(z .* b_i .* eps_i) / b_mix; % تقریب
        sig_mix = sum(z .* b_i .* sig_i) / b_mix; % تقریب

        % حل معادله حالت برای حجم مولی مخلوط
        f = @(V) P - R*T/(V - b_mix) + a_mix/((V + eps_mix*b_mix)*(V + sig_mix*b_mix));
        Vguess = R*T/P;
        V_m = fzero(f, Vguess);
        Z = P*V_m/(R*T);

        %% ضریب فوگاسیته برای متان
        i = 1; % شاخص متان
        Bi = b_i(i)*P/(R*T);
        B = b_mix*P/(R*T);
        A = a_mix*P/(R*T)^2;
        sum_aij = 0;
        for j=1:2
            sum_aij = sum_aij + z(j)*sqrt(a_i(i)*a_i(j));
        end
        delta = 2*sum_aij/a_mix - Bi/B;
        lnphi = Bi/B*(Z-1) - log((Z + sig_mix*B)/(Z - B)) ...
                + A/(B*(eps_mix - sig_mix))*delta*log((Z + eps_mix*B)/(Z + sig_mix*B));
        phi = exp(lnphi);

        %% حجم مولی جزئی از مشتق عددی فوگاسیته نسبت به فشار
        dP = 1e-3*P; % گام فشار
        f2 = @(V) (P+dP) - R*T/(V - b_mix) + a_mix/((V + eps_mix*b_mix)*(V + sig_mix*b_mix));
        V2 = fzero(f2, Vguess);
        Z2 = (P+dP)*V2/(R*T);
        lnphi2 = Bi/B*(Z2-1) - log((Z2 + sig_mix*B)/(Z2 - B)) ...
                + A/(B*(eps_mix - sig_mix))*delta*log((Z2 + eps_mix*B)/(Z2 + sig_mix*B));
        v_partial = R*T/(P) + R*T*(lnphi2 - lnphi)/dP;

        results = [results; T P Z phi v_partial];
    end
end

%% نمایش نتایج
fprintf('%6s %6s %6s %10s %10s\n','T(K)','P(bar)','Z','phi_CH4','V_m_partial');
for k=1:size(results,1)
    fprintf('%6.1f %6.1f %6.3f %10.3f %10.3f\n',results(k,1),results(k,2),results(k,3),results(k,4),results(k,5));
end
