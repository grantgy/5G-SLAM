function [x_post, P_poster]=kf_update(x_pre, P_pre, H, R, measurement)

    innov = measurement - H * x_pre;
    S = H * P_pre * H' + R;
    K = P_pre * H' / S;
    x_post = x_pre + K * innov;
    P_poster = (eye(size(P_pre)) - K * H) * P_pre;
    P_poster=(P_poster+P_poster')/2;
end