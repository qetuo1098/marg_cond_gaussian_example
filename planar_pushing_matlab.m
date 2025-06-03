% Author: Zi Cong Guo, zc.guo@mail.utoronto.ca

% Requires:
% Optimization Toolbox (for finding the constrained mean)

% Statistics and Machine Learning Toolbox (for generating noise)
% (can easily work around if not available - just replace 'mvnrnd')

% Image Processing Toolbox (only for plotting the probe)
% (can easily work around if not available - just replace 'viscircles')

%% 1. Initialize probe and box parameters
probe_radius = 0.1;
box_w = 0.2;
box_h = 1;
n_steps = 50;

probe_to_box_dist = probe_radius + box_w/2;
box_corners_box_frame = [-box_w/2, box_w/2, box_w/2, -box_w/2;
    box_h/2, box_h/2, -box_h/2, -box_h/2];
get_box_corners_f = @(box_pose) get_box_corners(box_pose, box_corners_box_frame);

%% 2. Generate groundtruth probe and box trajectory
% probe trajectory (noiseless)
probe_pos = [linspace(1, 3, n_steps);
    -linspace(0, 1.5, n_steps).^2];

% groundtruth box trajectory
box_coord_relative = [linspace(0, 0, n_steps);  % coord in sliding direction
    linspace(0, -pi/6, n_steps)];  % coord in rotating direction
box_coord_global = [probe_pos(1,:) + probe_to_box_dist;
    probe_pos(2,:) + box_coord_relative(1,:);
    box_coord_relative(2,:)];

% make the groundtruth box trajectory satisfy the constraints by bringing the box closer in the x direction
box_coord_global_adjust_x = zeros(1,n_steps);
opts = optimoptions(@fminunc,'Display','off');
for i = 1:n_steps
    box_corners = get_box_corners_f(box_coord_global(:,i));
    dist_to_probe = @(x) (p_poly_dist(probe_pos(1,i), probe_pos(2,i), box_corners(1,:) + x, box_corners(2,:)) - probe_radius)^2;
    box_coord_global_adjust_x(i) = fminunc(dist_to_probe, 0, opts);
end
box_coord_global(1,:) = box_coord_global(1,:) + box_coord_global_adjust_x;

%% 3. Generate box odom
Q = diag([0.03^2, 0.03^2, 0.01^2]);
u_noiseless = wrap_angle_in_state(box_coord_global(:,2:end) - box_coord_global(:,1:end-1));
u = wrap_angle_in_state(u_noiseless + mvnrnd([0 0 0], Q, n_steps-1)');  % add noise to the odom

% assume no measurements for this script, so we're just dead reckoning.
% one can implement the contact measurements in the InCOpt scenario if desired. The results are visually similar.

%% 4. Optimize for the constrained mean
objFn = @(x) compute_cost(x, box_coord_global(:,1), u, Q);
constr = @(x) box_touch_probe_constraint_wrapper(x, probe_pos(:,2:end), probe_radius, get_box_corners_f);
constr_eq = @(x) box_touch_probe_constraint(x, probe_pos(:,2:end), probe_radius, get_box_corners_f);
x_init = box_coord_global(:,2:end);
opts = optimoptions(@fmincon,'MaxFunEvals',5e4);

% constrained localization - find the optimal constrained trajectory
[x, fval, exitflag, output, lambda, grad, hessian] = fmincon(objFn, x_init, [], [], [], [], [], [], constr, opts);
% Doing this brute force with the MATLAB solver is slow, but one can easily speed this up by
% exploiting sparsity by swapping this brute-force solver with, e.g., InCOpt.

%% 5. Compute the unconstrained covariance
% Compute the unconstrained covariance through the Gauss-Newton.
% We approximate of the cost function's Hessian: F = H^T W_inv H, where H is the Jacobian.
% For details, see [T. D. Barfoot, State Estimation for Robotics.
% Cambridge, U.K.: Cambridge Univ. Press, 2017.]
H = zeros((n_steps-1)*3, (n_steps-1)*3);
for i = 1:n_steps-1
    H(i*3-2:i*3, i*3-2:i*3) = -eye(3);
    if i < n_steps-1
        H(i*3+1:i*3+3, i*3-2:i*3) = eye(3);
    end
end
WCell = repmat({Q}, 1, n_steps-1);
unconstrained_info = H' / blkdiag(WCell{:}) * H;
unconstrained_cov = inv(unconstrained_info);
% by inspecting the diagonal blocks of unconstrained_cov, we can see that
% the covariance is  a perfect sphere that increases with every timestep. This makes sense, 
% since we're simply adding an isotropic Gaussian noise at every step and we're just dead reckoning

%% 6. compute the constrained covariance
% note the simplicity: just 3 lines of code!
S = numdiff(constr_eq, x);  % compute linearized constraints
N = null(S');
constrained_cov = N / (N' * unconstrained_info * N) * N';  % (II.4) in paper
% If desired, one can exploit sparsity to speed up the computation of null(S') and inv(N' * unconstrained_info * N).
% For an example, see Appendix E in [Z. C. Guo, F. Dümbgen, J. R. Forbes, and T. D. Barfoot, "Data-driven
% batch localization and SLAM using Koopman linearization,” IEEE T-RO, vol. 40, pp. 3964–3983, 2024.]

%% not done here: map back to the manifold
% This is mapping is not always necessary depending why we need the
% covariance. It is necessary for evaluating consistency, or else the Mahalanobis
% distance is not computable. For this scenario, we can realize that the structure is
% actually a revolute joint followed by a prismatic joint. See footnote 7
% in paper, and also kinematic_chain.png in this github repo.

%% 7. plot
% Here we plot the pose covariance of the full batch solution at every timestep.

figure
plot(box_coord_global(1,:), box_coord_global(2,:), 'Color', 'black', 'DisplayName', 'Groundtruth');
axis equal
hold on
legend
plot(x(1,:), x(2,:), 'Color', 'Blue', 'DisplayName', 'Mean estimate')
for i = 1:n_steps-1
    box_corners_gt = get_box_corners_f(box_coord_global(:,i+1));
    plt_box_gt = plot(polyshape(box_corners_gt(1,:), box_corners_gt(2,:)), 'FaceColor','Black', 'DisplayName', 'Groundtruth');
    box_corners_est = get_box_corners_f(x(:,i));
    plt_box_est = plot(polyshape(box_corners_est(1,:), box_corners_est(2,:)), 'FaceColor','Blue', 'DisplayName', 'Mean estimate');
    plt_probe = viscircles(probe_pos(:,i+1)', probe_radius, 'Color', 'Blue');
    plt_unc_cov = plotcov(x(1:2,i), unconstrained_cov(i*3-2:i*3-1, i*3-2:i*3-1), 3, 'Color', 'Red', 'DisplayName', 'Unconstrained Cov.');
    plt_con_cov = plotcov(x(1:2,i), constrained_cov(i*3-2:i*3-1, i*3-2:i*3-1), 3, 'Color', 'Green', 'DisplayName', 'Constrained Cov.');
    pause(0.2)
    if i < n_steps-1
        delete(plt_box_gt)
        delete(plt_box_est)
        delete(plt_probe)
        delete(plt_unc_cov)
        delete(plt_con_cov)
    end
end
% The ellipse corresponding to the constrained covariance is very thin, but it still has some
% width - it's not simply a line. This is because the box can rotate around
% the probe, which means the box pose can move slightly off the sliding line.

% This makes sense dimensions-wise: since the constraint (must be in contact) is 1D,
% the box has 3-1=2 degrees of freedom (sliding and rotating). As a result, the covariance of its 
% pose (x,y,theta) is rank 2, which is degenerate. But when plotting its x-y covariance ellipse above,
% we're marginalizing that covariance to only (x,y). This marginalized covariance is still rank 2,
% meaning it is not degenerate, so the ellipse is not a line but has some width.

% Small note: this is slightly different than in the ICRA presentation on Youtube (and in the paper),
% which instead plots the end-pose covariance at each timestep given the trajectory so far, since it
% was doing incremental localization. This means that at the time of plotting, the poses at the later
% timesteps were not available.
% Visually, though, there's not much difference between this script and the one on YouTube (and in the paper).
%% Helper functions
function x = wrap_angle_in_state(x)
x(3,:) =  wrapToPi(x(3,:));
end

function box_corners = get_box_corners(box_pose, box_corners_box_frame)
rot_mat = rotz( rad2deg(box_pose(3)));
box_corners = rot_mat(1:2,1:2) * box_corners_box_frame + box_pose(1:2);
end

function residual = compute_residual(x,u)
predicted_state = wrap_angle_in_state(x(:,1:end-1) + u);
residual = wrap_angle_in_state(x(:,2:end) - predicted_state);
end

function cost = compute_cost(x, x0, u, Q)
% cost is just the motion prior
residual = compute_residual([x0 x], u);
cost = 0;
for i = 1:size(residual, 2)
    cost = cost + residual(:,i)' / Q * residual(:,i);
end
end

function ceq = box_touch_probe_constraint(x, probe_pos, probe_radius, get_box_corners_f)
ceq = zeros(1, size(x,2));
for i = 1:size(x,2)
    box_corners = get_box_corners_f(x(:,i));
    ceq(i) = (p_poly_dist(probe_pos(1,i), probe_pos(2,i), box_corners(1,:), box_corners(2,:)) - probe_radius)^2;
end

end
function [c, ceq] = box_touch_probe_constraint_wrapper(x, probe_pos, probe_radius, get_box_corners_f)
c = [];
ceq = box_touch_probe_constraint(x, probe_pos, probe_radius, get_box_corners_f);
end

function S = numdiff(fn, x)
% very simple numerical differentiation using central difference
eps = 1e-7;
n_i = size(x, 1);
n_j = size(x, 2);
val = reshape(fn(x), [], 1);
S = zeros(size(val,1), n_i*n_j);
for i = 1:n_i*n_j
    h = zeros(n_i, n_j);
    h(i) = eps;
    S(:, i) = reshape((fn(x + h) - fn(x - h)) / (2*eps), [], 1);
end
S = S';
end
