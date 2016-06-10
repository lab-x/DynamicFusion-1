
#include <iostream>
#include <cvd/image_io.h>
#include <cvd/image_ref.h>
#include <Eigen/CholmodSupport>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>



using namespace std;

int k_n = 3;
float eps = 0.25;
int height = 240;
int width = 320;
float er = 0.0001;
int _frame = 0;

struct n_i {
	Eigen::Vector3d v;
	float w;
	Eigen::Matrix4d se3;
};

std::vector<n_i> n_warp;

void findKnearestPointIndex(int* out, std::vector<n_i>& n_warp, Eigen::Vector3d x, int index) {
	int n_size = n_warp.size();
	int checked = 0;
	unsigned int k;
	for (k = 0; k < k_n; ++k) {
		if (k >= n_size) {
			out[k] = -1;
			continue;
		}
		int min = 0;
		float min_d = std::numeric_limits<float>::max();

		for(uint i = 0; i < n_size; ++i) {
			if (index == i) {
				continue;
			}
			float d = (n_warp[i].v - x).norm();
			if (d < min_d && (checked >> i) % 2 == 0) {
				min = i;
				min_d = d;
			}
		}
		out[k] = min;
		checked |= 1 << min;
	}
}

void getWarpMatrix(Eigen::Matrix4d* out, Eigen::Vector3d x) {
	if (n_warp.size() == 0) {
		*out = Eigen::Matrix4d::Identity();
		return;
	}

	//Eigen::VectorXd bq;
	Eigen::Matrix4d bm = Eigen::Matrix4d::Zero();
	// Get k nearest points
	int* kNear = new int[k_n];
	findKnearestPointIndex(kNear, n_warp, x, -1);
	// blend them
	double weight_sum = 0;
	unsigned int k;
	for (k = 0; k < k_n; ++k) {
      /*
		// For BDQ
		Eigen::Quaternion qr(n_warp[k].se3.block(0,0,3,3));
		Eigen::Vector4d t;
		t << 1 << n_warp[k].v ;
		Eigen::Quaternion qt(t);
		Eigen::VectorXd dq;
		dq << qr << 0.5*qt * qr;
		bq += n_warp[k].w * dq;
	}
	bq /= bq.norm();
	Eigen::Quaternion qr(bq.head(4));
	Eigen::Matrix3d r = qr.matrix();
	bq.tail(4)
	*/
		if (kNear[k] < 0) {
			break;
		}
		bm += n_warp[kNear[k]].w * n_warp[kNear[k]].se3;
		weight_sum += n_warp[kNear[k]].w;
	}
	*out = bm/weight_sum;
}

void vertex2normalKernel(Eigen::Vector3d * out, const Eigen::Vector3d * in) {
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			const Eigen::Vector3d left = in[max(x - 1, 0) + width * y];
			const Eigen::Vector3d right = in[min(x + 1, width - 1) + width * y];
			const Eigen::Vector3d up = in[x + width * max(int(y) - 1, 0)];
			const Eigen::Vector3d down = in[x + width * min(y + 1, ((int) height) - 1)];

			if (left(2) == 0 || right(2) == 0 || up(2) == 0 || down(2) == 0) {
				out[x + y * width] = Eigen::Vector3d(0,0,0);
				continue;
			}
			const Eigen::Vector3d dxv = right - left;
			const Eigen::Vector3d dyv = down - up;

			Eigen::Vector3d normal = dyv.cross(dxv);
			normal.normalize();
			out[x + y * width] = normal; // switched dx and dy to get factor -1
		}
	}
}

void so3Matrix(Eigen::Matrix3d* rm, Eigen::Vector3d v) {
	*rm << 1, -v[2], v[1],
	v[2], 1, -v[0],
	-v[1], v[0], 1;
}

bool non_rigid_track(Eigen::Vector3d* vertex, Eigen::Vector3d* normal,
					Eigen::Matrix4d cameraMatrix,
					Eigen::Vector3d* inputVertex,
					Eigen::Matrix4d pose)
{

	if (n_warp.size() == 0) {
		return true;
	}


	// Reg term
	int n_dp = n_warp.size();

	Eigen::SparseMatrix<double> A_reg(n_dp * k_n * 3, n_dp * 6);
	Eigen::VectorXd b_reg(n_dp * k_n * 3);
	A_reg.setZero();

	Eigen::Matrix4d rigid_t = pose;

	//std::cout << "rigid_t" << rigid_t << std::endl;
	unsigned int i;

	for (i = 0; i < n_dp; ++i) {
		int* kNear = new int[k_n];

		/// For each deformation node, we are trying to search for the
		/// vertices near to it.

		findKnearestPointIndex(kNear, n_warp, n_warp[i].v, i);

		for (int k = 0; k < k_n; k++)
		{
			if (kNear[k] < 0) {
				return false;
			}

			int j = kNear[k];

			/// find the max of the weight of two vertices
			double alpha = n_warp[j].w > n_warp[i].w ? n_warp[j].w : n_warp[i].w;
			double con_a = pow(alpha, 0.5);

			Eigen::Matrix4d non_rigid_t_j;
			getWarpMatrix(&non_rigid_t_j, n_warp[j].v);
			Eigen::Matrix4d non_rigid_t_i;
			getWarpMatrix(&non_rigid_t_i, n_warp[i].v);

			Eigen::Vector3d ridgv =  non_rigid_t_i.block(0,0,3,3) * n_warp[j].v;

			//std::cout << "rigid_t.block(0,0,3,3) " << k << "  " << rigid_t.block(0,0,3,3)  << std::endl;

			A_reg.coeffRef((i * k_n + k)*3, i * 6) = 0;
			A_reg.coeffRef((i * k_n + k)*3 + 1, i * 6)      = con_a * ridgv[2];
			A_reg.coeffRef((i * k_n + k)*3 + 2, i * 6)      = -con_a * ridgv[1];

			A_reg.coeffRef((i * k_n + k)*3, i * 6 + 1)      = -con_a * ridgv[2];
			A_reg.coeffRef((i * k_n + k)*3 + 1, i * 6 + 1)  = 0;
			A_reg.coeffRef((i * k_n + k)*3 + 2, i * 6 + 1)  = con_a * ridgv[0];

			A_reg.coeffRef((i * k_n + k)*3, i * 6 + 2)      = con_a * ridgv[1];
			A_reg.coeffRef((i * k_n + k)*3 + 1, i * 6 + 2)  = -con_a * ridgv[0];
			A_reg.coeffRef((i * k_n + k)*3 + 2, i*6 + 2)    = 0;

			A_reg.coeffRef((i * k_n + k)*3, i * 6 + 3)      = con_a;
			A_reg.coeffRef((i * k_n + k)*3 + 1, i * 6 + 4)  = con_a;
			A_reg.coeffRef((i * k_n + k)*3 + 2, i * 6 + 5)  = con_a;

			Eigen::Vector3d rjdgv = non_rigid_t_j.block(0,0,3,3) * n_warp[j].v;

			A_reg.coeffRef((i * k_n + k)*3, j * 6)          = 0;
			A_reg.coeffRef((i * k_n + k)*3 + 1, j * 6)      = -con_a * rjdgv[2];
			A_reg.coeffRef((i * k_n + k)*3 + 2, j * 6)      = con_a * rjdgv[1];

			A_reg.coeffRef((i * k_n + k)*3, j * 6 + 1)      = con_a * rjdgv[2];
			A_reg.coeffRef((i * k_n + k)*3 + 1, j * 6 + 1)  = 0;
			A_reg.coeffRef((i * k_n + k)*3 + 2, j*6 + 1)    = -con_a * rjdgv[0];

			A_reg.coeffRef((i * k_n + k)*3, j * 6 + 2)      = -con_a * rjdgv[1];
			A_reg.coeffRef((i * k_n + k)*3 + 1, j * 6 + 2)  = con_a * rjdgv[0];
			A_reg.coeffRef((i * k_n + k)*3 + 2, j * 6 + 2)  = 0;

			A_reg.coeffRef((i * k_n + k)*3, j * 6 + 3)      = -con_a;
			A_reg.coeffRef((i * k_n + k)*3 + 1, j * 6 + 4)  = -con_a;
			A_reg.coeffRef((i * k_n + k)*3 + 2, j * 6 + 5)  = -con_a;

			Eigen::Vector3d ti =  non_rigid_t_i.block(0,3,3,1);
			Eigen::Vector3d tj =  non_rigid_t_j.block(0,3,3,1);
			Eigen::Vector3d t =  con_a * (ridgv - rjdgv + ti - tj);

			b_reg((i * k_n + k)*3)     = t[0];
			b_reg((i * k_n + k)*3 + 1) = t[1];
			b_reg((i * k_n + k)*3 + 2) = t[2];
		}
	}
	// Data term
	Eigen::SparseMatrix<double> A_data(width * height, n_dp * 6);
	Eigen::VectorXd b_data(width * height);
	A_data.setZero();
	int row = 0;
	int pixelx, pixely;
	for (pixely = 0; pixely < height; pixely++) {
		for (pixelx = 0; pixelx < width; pixelx++) {
			Eigen::Vector3d v_u = vertex[pixelx + pixely * width];
			Eigen::Vector3d r_v_u = rigid_t.block(0,0,3,3) * v_u + rigid_t.block(0,3,3,1);
			Eigen::Vector3d n_u = normal[pixelx + pixely * width];

			if (n_u == Eigen::Vector3d(0,0,0)) {
				continue;
			}

			Eigen::Vector3d r_n_u = rigid_t.block(0,0,3,3) * n_u;

			Eigen::Matrix4d non_rigid_t;
			getWarpMatrix(&non_rigid_t, v_u);

			Eigen::Matrix4d eigenCameraMatrix = cameraMatrix;

			Eigen::Vector3d v_hat_u = non_rigid_t.block(0,0,3,3) * r_v_u + non_rigid_t.block(0,3,3,1); // Non_rigid warp is initialised to be Identity
			Eigen::Vector3d projected_u_hat = eigenCameraMatrix.block(0,0,3,3) * v_hat_u + eigenCameraMatrix.block(0,3,3,1);
			//Eigen::Vector3d projected_u_hat = eigenCameraMatrix.block(0,0,3,3) * (w.block(0,0,3,3) * v_u + w.block(0,3,3,1))
			//		+ eigenCameraMatrix.block(0,3,3,1);
			Eigen::Vector2i pixel_u;

			pixel_u << int(projected_u_hat[0] / projected_u_hat[2] + 0.5),
					int(projected_u_hat[1] / projected_u_hat[2] + 0.5);
			if (pixel_u[0] < 0 || pixel_u[0] > width - 1
					|| pixel_u[1] < 0 || pixel_u[1] > height - 1) {
				continue;
			}

			Eigen::Vector3d vl = inputVertex[pixel_u[0] + pixel_u[1] * width];

			if (vl.norm() < er) {
				continue;
			}


			int* kNear = new int[k_n];
			findKnearestPointIndex(kNear, n_warp, v_u, -1);

			Eigen::Vector3d C_n(0, 0, 0);
			Eigen::Vector3d C_v(0, 0, 0);
			double eta = 0;

			for (int k = 0; k < k_n; k++) {
				n_i d_p = n_warp[kNear[k]];
				C_n += d_p.w * d_p.se3.block(0,0,3,3) * r_n_u;
				C_v += d_p.w * d_p.se3.block(0,0,3,3) * r_v_u + d_p.se3.block(0,3,3,1);
				eta += d_p.w;
			}

			Eigen::Vector3d C_n_t = C_n;
			C_n /= eta;
			C_n.normalize();
			C_v /= eta;

			Eigen::Vector3d D_v = C_v - vl;

			for (int k = 0; k < k_n; k++) {
				n_i d_p = n_warp[kNear[k]];

				Eigen::Vector3d rn = d_p.se3.block(0,0,3,3) * r_n_u;
				rn.normalize();
				Eigen::Matrix3d rnm;
				so3Matrix(&rnm, rn);

				Eigen::Vector3d rv = d_p.se3.block(0,0,3,3) * r_v_u + d_p.se3.block(0,3,3,1);
				Eigen::Matrix3d rvm;
				so3Matrix(&rvm, rv);

				Eigen::RowVector3d c_r = d_p.w * (D_v.transpose() * (rnm) +
						C_n.transpose() * (rvm)) / eta;

				Eigen::RowVector3d c_t = d_p.w / eta * C_n.transpose();

				int j = kNear[k];
				A_data.coeffRef(row, j * 6)     = c_r[0];
				A_data.coeffRef(row, j * 6 + 1) = c_r[1];
				A_data.coeffRef(row, j * 6 + 2) = c_r[2];

				A_data.coeffRef(row, j * 6 + 3) = c_t[0];
				A_data.coeffRef(row, j * 6 + 4) = c_t[1];
				A_data.coeffRef(row, j * 6 + 5) = c_t[2];
			}

			b_data(row) = C_n.transpose() * D_v;
			if (isnan(C_n.transpose() * D_v)) {
				std::cout << "C_n_t" << C_n_t << std::endl;
				std::cout << "eta" << eta << std::endl;
				std::cout << "r_n_u" << r_n_u << std::endl;
				//std::cout << "D_v" << D_v << std::endl;

			}
			row ++;
		}
	}

	A_data.conservativeResize(row,n_dp * 6);
	b_data.conservativeResize(row);

	float w_data = 0.5;
	float w_reg = 100;
	float w_lm = 1E-12;

	float reg = 0;

	for (int i = 0; i < n_dp; i ++) {
		int* kNear = new int[k_n];
		findKnearestPointIndex(kNear, n_warp, n_warp[i].v, i);
		for (int k = 0; k < k_n; k++) {
			int j = kNear[k];
			double alpha = n_warp[j].w > n_warp[i].w ? n_warp[j].w : n_warp[i].w;
			Eigen::Matrix4d Ti;
			getWarpMatrix(&Ti, n_warp[i].v);
			Eigen::Matrix4d Tj;
			getWarpMatrix(&Tj, n_warp[j].v);
			Eigen::Vector3d tidgv =  Ti.block(0,0,3,3) * n_warp[j].v + Ti.block(0,3,3,1);
			Eigen::Vector3d tjdgv =  Tj.block(0,0,3,3) * n_warp[j].v + Tj.block(0,3,3,1);
			Eigen::Vector3d v = tidgv - tjdgv;
			reg += 0.5 * alpha * ((v.transpose() * v)[0]);
		}
	}

	std::cout << "reg " << reg;

	float data = 0;

	for (int pixely = 0; pixely < height; pixely++) {
		for (int pixelx = 0; pixelx < width; pixelx++) {
			Eigen::Vector3d v_u  = vertex[pixelx + pixely * width];
			Eigen::Vector3d r_v_u = rigid_t.block(0,0,3,3) * v_u + rigid_t.block(0,3,3,1);
			Eigen::Vector3d n_u = normal[pixelx + pixely * width];
			Eigen::Vector3d r_n_u = rigid_t.block(0,0,3,3) * n_u;

			Eigen::Matrix4d non_rigid_t;
			getWarpMatrix(&non_rigid_t, v_u);


			Eigen::Matrix4d eigenCameraMatrix = cameraMatrix;
			Eigen::Vector3d v_hat_u = non_rigid_t.block(0,0,3,3) * r_v_u + non_rigid_t.block(0,3,3,1);
			Eigen::Vector3d projected_u_hat = eigenCameraMatrix.block(0,0,3,3) * v_hat_u + eigenCameraMatrix.block(0,3,3,1);
			Eigen::Vector3d n_hat_u = non_rigid_t.block(0,0,3,3) * r_n_u;

			Eigen::Vector2i pixel_u;

			pixel_u << int(projected_u_hat[0] / projected_u_hat[2] + 0.5),
					int(projected_u_hat[1] / projected_u_hat[2] + 0.5);
			if (pixel_u[0] < 0 || pixel_u[0] > width - 1
					|| pixel_u[1] < 0 || pixel_u[1] > height - 1) {
				continue;
			}
			Eigen::Vector3d vl = inputVertex[pixel_u[0] + pixel_u[1] * width];
			data += 0.5 * pow((n_hat_u.transpose() * (v_hat_u - vl))[0], 2);
		}
	}

	std::cout << "data " << data;

	float energy = reg + w_reg * 2.0 * data;
	std::cout << " This is energy ---------------------->" << energy << std::endl;


	Eigen::SparseMatrix<double>A_lm(n_dp * 6, n_dp * 6);

	for(int i = 0; i < n_dp * 6; i++)
	{
	    A_lm.insert(i,i) = 1;
	}

	Eigen::CholmodSimplicialLDLT< Eigen::SparseMatrix<double> >solver;
	solver.compute(w_reg * A_reg.transpose()*A_reg + w_data * A_data.transpose()*A_data + w_lm * A_lm);

	if(solver.info() != Eigen::Success)
	{
	    // decomposition failed
	    std::cout<<"Decomposition failed" << std::endl;
	    return false;
	}

	//std::cout << "solved, - n_dp:" << n_dp << std::endl;
	Eigen::VectorXd x_update;
	Eigen::VectorXd Axb = -1.0f * (w_reg * A_reg.transpose() * b_reg + w_data * A_data.transpose() * b_data);
	x_update = solver.solve(Axb);

    std::cout<<"norm of x_update = " << x_update.norm()/ (n_dp) << std::endl;

	//TODO Set a timer here
	//std::cout<<"solver time elapsed = " << solver_timer.elapsed() << std::endl;

	// Updates
	for(int i = 0 ; i < n_dp;i++)
	{
	    // Get the rotation update
	    Eigen::Vector3d so3_i(x_update(6*i+0), x_update(6*i+1), x_update(6*i+2));
	    //std::cout<<"so3_i"<< so3_i << std::endl;
	    // Get the translation update
	    Eigen::Vector3d t_i(x_update(6*i+3), x_update(6*i+4), x_update(6*i+5));

	    Eigen::Matrix3d r_i;
  	    so3Matrix(&r_i, so3_i);

	    Eigen::Matrix3d o_r = n_warp[i].se3.block(0,0,3,3);
	    Eigen::Vector3d o_t = n_warp[i].se3.block(0,3,3,1);

	    Eigen::Matrix3d n_r = r_i * o_r;
	    Eigen::Vector3d n_t = o_t + t_i;

	    Eigen::Matrix4d n_w;
	    n_w << n_r, n_t,
	    	   0,0,0,1;
	    n_warp[i].se3 = n_w;
	    //std::cout<<"pose("<<i<<") = "<< n_w << std::endl;
	}
	return true;
}

bool test_non_rigid_track(Eigen::Vector3d* vertex, Eigen::Matrix4d cameraMatrix, Eigen::Vector3d* inputVertex, Eigen::Matrix4d pose) {
	Eigen::Matrix4d rigid_t = pose;
	double avg_diff_norm = 0.0;
	CVD::Image<u_int16_t>depthmap_n(CVD::ImageRef(width, height));
	int unacceptable_n = 0;
	int pixelx, pixely;
	for (pixely = 0; pixely < height; pixely++) {
		for (pixelx = 0; pixelx < width; pixelx++) {
			depthmap_n[CVD::ImageRef(pixelx,pixely)] = (u_int16_t)0;

			Eigen::Vector3d v_u_k = vertex[pixelx + pixely * width];

			Eigen::Matrix4d w; // = Eigen::Matrix4d::Identity();
			getWarpMatrix(&w, v_u_k);
			v_u_k = pose.block(0,0,3,3) * v_u_k + pose.block(0,3,3,1);

			v_u_k = w.block(0,0,3,3) * v_u_k + w.block(0,3,3,1);
			Eigen::Matrix4d eigenCameraMatrix = cameraMatrix;
			Eigen::Vector3d projected_u_hat = eigenCameraMatrix.block(0,0,3,3) * v_u_k + eigenCameraMatrix.block(0,3,3,1);
			Eigen::Vector2i pixel_u;
			pixel_u << int(projected_u_hat[0] / projected_u_hat[2] + 0.5),
					int(projected_u_hat[1] / projected_u_hat[2] + 0.5);
			if (pixel_u[0] < 0 || pixel_u[0] > width - 1
					|| pixel_u[1] < 0 || pixel_u[1] > height - 1) {
				continue;
			}

			Eigen::Vector3d vl = inputVertex[pixel_u[0] + pixel_u[1] * width];
			Eigen::Vector3d diff = v_u_k - vl;
			double norm = diff.norm();
			depthmap_n[CVD::ImageRef(pixelx,pixely)] = norm * 5000 ;

			avg_diff_norm += norm;
		}
	}

	char fileName[200];
	std::sprintf(fileName,"test_diff_%05d.png",_frame);
	CVD::img_save(depthmap_n, fileName);
	std::cout << "average test diff: " << avg_diff_norm/(height * width) << std::endl;
}

int main() {
	float depthMap[width * height];

	Eigen::Matrix4d K;
	K << 481.20 /2, 0.00, 319.50 /2, 0,
		 0.00, -480.00 / 2,	239.50 / 2, 0,
		 0.00, 0.00, 1.00, 0,
		 0, 0, 0, 1;

	Eigen::Vector3d refvertex[height * width];
	Eigen::Vector3d invertex[height * width];
	Eigen::Vector3d normal[height * width];
	Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();


	for (int i = 0; i < 50; i++) {
		char fileNamepng[200];
		sprintf(fileNamepng,"%s/test_%05d.png","../data",_frame);
		CVD::Image<u_int16_t>DEPTHPNG;
		CVD::img_load(DEPTHPNG, fileNamepng);

		for (int v = 0; v < height ; v++)
		{
			for (int u = 0; u < width; u++)
			{
				depthMap[u + v * width] = 0;

			}
		}

		std::cout << "reset depth map" << std::endl;

		for (int v = 0; v < height; v++)
		{
			for (int u = 0; u < width; u++)
			{
				bool inrange = DEPTHPNG[CVD::ImageRef(u,v)] < 8000 && u > 21 && u < 200;
				//bool inrange = true;
				depthMap[u + v * width] = inrange ? ((float)DEPTHPNG[CVD::ImageRef(u,v)])/5000.0f : 0;

			}
		}

		std::cout << "get depth map" << std::endl;

		n_warp.clear();

		for (int v = 0; v < height ; v++)
		{
			for (int u = 0; u < width; u++)
			{
				Eigen::Matrix4d invK = K.inverse();
				Eigen::Vector3d v_p(v,u,1);
				Eigen::Vector3d v_c = depthMap[u + v * width] * (invK.block(0,0,3,3) * v_p + invK.block(0,3,3,1));
				invertex[u + v * width] = v_c;

				if (abs(v_c(2) - 0.0) <= er  && abs(v_c(1) - 0.0) <= er && abs(v_c(0) - 0.0) <= er) {
					continue;
				}
				int* kNear = new int[k_n];
				findKnearestPointIndex(kNear, n_warp, v_c, -1);
				double min_dis = std::numeric_limits<double>::max();
				for (int k = 0; k < k_n; k++) {
					if (kNear[k] < 0) {
						continue;
					}

					double dis = (n_warp[kNear[k]].v - v_c).norm() / n_warp[kNear[k]].w;
					min_dis = dis < min_dis ? dis : min_dis;
				}
				if (min_dis >=1) {
					n_i dgnew;
					dgnew.v = v_c;
					Eigen::Matrix4d wm;
					getWarpMatrix(&wm, v_c);
					dgnew.se3 = wm;
					dgnew.w = eps;
					n_warp.push_back(dgnew);
					std::cout << "added "<< "  " << dgnew.v  << std::endl;
					//std::cout << "added se3"<< "  " << dgnew.se3  << std::endl;
				}

			}
		}


		if (i > 0)
		{
			bool non_rigid_return = false;
			for (int n = 0; n < 6; n++) {
				non_rigid_return = non_rigid_track(refvertex, normal,K, invertex, pose);
				test_non_rigid_track(refvertex, K, invertex, pose);
			}
			std::cout << "non_rigid tracking end" << std::endl;
		}

		vertex2normalKernel(normal, invertex);

		for (int v = 0; v < height; v++)
		{
			for (int u = 0; u < width; u++)
			{
				refvertex[u + v * width] = invertex[u + v * width];
			}
		}

		_frame++;
	}

	return 0;
}
