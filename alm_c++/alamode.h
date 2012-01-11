namespace ALM_NS {
	class ALM {
	public:
		class Memory *memory; // Memory Allocation
		class Error *error; // Error handling
		class Timer *timer; // CPU timer

		ALM(int, char **);
		~ALM();
	};
}