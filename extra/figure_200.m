figure();
post1 = load('post_200_1.mat');
post2 = load('post_200_2.mat');
post3 = load('post_200_3.mat');
load('post_vb.mat')
sam1 = load('samples_200_1.mat');
sam2 = load('samples_200_2.mat');
sam3 = load('samples_200_3.mat');
N = 21;
angle = linspace(0,pi*(N-1)/N,N);
fix_item = fix(N/2)+1;
% subplot(2,2,[1,3])
plot(angle,exp(post1.p_final-logsumexp(post1.p_final)),'-o','LineWidth', 2);
hold on;
plot(angle,exp(post2.p_final-logsumexp(post2.p_final)),'-o','LineWidth', 2);
hold on;
plot(angle,exp(post3.p_final-logsumexp(post3.p_final)),'-o','LineWidth', 2);


figure();
% subplot(2,2,2)
plot(angle,sum(sam1.x_samples,2),'-*','LineWidth',2);
hold on;
plot(angle,sum(sam2.x_samples,2),'-*','LineWidth',2);
hold on;
plot(angle,sum(sam3.x_samples,2),'-*','LineWidth',2);

