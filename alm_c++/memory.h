namespace ALM_NS {
	class Memory {
	public:
		Memory();
		~Memory();

		double **create_d2_array(int, int);
		void destroy_d2_array();
		double **create_i2_array(int, int);
		void destroy_i2_array();

	};
}