pragma solidity ^0.8.0;

contract MultiSigWallet {
    uint public numConfirmationsRequired;
    address[] public owners;
    mapping(address => bool) public isOwner;
    uint public transactionCount;
    mapping(uint => Transaction) public transactions;
    mapping(uint => mapping(address => bool)) public isConfirmed;
    uint256 public constant MAX_OWNER_COUNT = 50;
    uint256 public lastActionTimestamp;
    uint256 public immutable CONFIRMATION_TIMEOUT;
    uint256 public dailyLimit;
    uint256 public spentToday;
    uint256 public lastDay;
    bool public paused;
    
    struct Transaction {
        address to;
        uint value;
        bytes data;
        bool executed;
        uint numConfirmations;
        uint256 submissionTime;
        address submitter;
        bool canceled;
        string description;
        uint256 nonce;
    }
    
    struct DailyOperation {
        uint256 value;
        uint256 timestamp;
        address initiator;
        bool processed;
    }
    
    mapping(uint256 => DailyOperation) public dailyOperations;
    uint256 public dailyOperationCount;
    mapping(address => uint256) public lastConfirmationTime;
    mapping(address => uint256) public ownerActivityCount;
    mapping(bytes32 => bool) public usedNonces;
    
    event Deposit(address indexed sender, uint value, uint256 timestamp);
    event Submission(uint indexed transactionId, address indexed submitter);
    event Confirmation(address indexed sender, uint indexed transactionId, uint256 timestamp);
    event Execution(uint indexed transactionId, uint256 executionTime);
    event ExecutionFailure(uint indexed transactionId, string reason);
    event OwnerAddition(address indexed owner);
    event OwnerRemoval(address indexed owner);
    event DailyLimitChange(uint256 dailyLimit);
    event TransactionCanceled(uint indexed transactionId);
    event ContractPaused(address indexed initiator);
    event ContractUnpaused(address indexed initiator);
    event ConfirmationTimeoutChanged(uint256 newTimeout);
    
    modifier notNull(address _address) {
        require(_address != address(0), "Invalid null address");
        _;
    }
    
    modifier onlyOwner() {
        require(isOwner[msg.sender], "Not an owner");
        _;
    }
    
    modifier ownerDoesNotExist(address owner) {
        require(!isOwner[owner], "Owner already exists");
        _;
    }
    
    modifier ownerExists(address owner) {
        require(isOwner[owner], "Owner does not exist");
        _;
    }
    
    modifier transactionExists(uint _transactionId) {
        require(transactions[_transactionId].to != address(0), "Transaction does not exist");
        _;
    }
    
    modifier confirmed(uint _transactionId) {
        require(isConfirmed[_transactionId][msg.sender], "Transaction not confirmed");
        _;
    }
    
    modifier notConfirmed(uint _transactionId) {
        require(!isConfirmed[_transactionId][msg.sender], "Transaction already confirmed");
        _;
    }
    
    modifier notExecuted(uint _transactionId) {
        require(!transactions[_transactionId].executed, "Transaction already executed");
        _;
    }
    
    modifier notCanceled(uint _transactionId) {
        require(!transactions[_transactionId].canceled, "Transaction was canceled");
        _;
    }
    
    modifier validRequirement(uint ownerCount, uint _required) {
        require(ownerCount <= MAX_OWNER_COUNT && _required <= ownerCount && _required > 0);
        _;
    }
    
    modifier whenNotPaused() {
        require(!paused, "Contract is paused");
        _;
    }
    
    constructor(address[] memory _owners, uint _numConfirmationsRequired) 
        validRequirement(_owners.length, _numConfirmationsRequired)
    {
        require(_owners.length > 0, "Owners required");
        require(
            _numConfirmationsRequired > 0 && _numConfirmationsRequired <= _owners.length,
            "Invalid number of required confirmations"
        );
        
        for (uint i = 0; i < _owners.length; i++) {
            address owner = _owners[i];
            require(owner != address(0), "Invalid owner");
            require(!isOwner[owner], "Duplicate owner");
            
            isOwner[owner] = true;
            owners.push(owner);
            ownerActivityCount[owner] = 0;
            lastConfirmationTime[owner] = block.timestamp;
        }
        
        numConfirmationsRequired = _numConfirmationsRequired;
        CONFIRMATION_TIMEOUT = 24 hours;
        lastActionTimestamp = block.timestamp;
        dailyLimit = 100 ether;
        spentToday = 0;
        lastDay = block.timestamp / 1 days;
        paused = false;
    }
    
    receive() external payable {
        if (msg.value > 0) {
            emit Deposit(msg.sender, msg.value, block.timestamp);
        }
    }
    
    function submitTransaction(
        address _to,
        uint _value,
        bytes memory _data,
        string memory _description
    ) 
        public
        onlyOwner
        whenNotPaused
        notNull(_to)
        returns (uint transactionId)
    {
        require(_value <= address(this).balance, "Insufficient balance");
        
        transactionId = transactionCount++;
        Transaction storage transaction = transactions[transactionId];
        transaction.to = _to;
        transaction.value = _value;
        transaction.data = _data;
        transaction.executed = false;
        transaction.numConfirmations = 0;
        transaction.submissionTime = block.timestamp;
        transaction.submitter = msg.sender;
        transaction.canceled = false;
        transaction.description = _description;
        transaction.nonce = uint256(keccak256(abi.encodePacked(block.timestamp, msg.sender, transactionId)));
        
        emit Submission(transactionId, msg.sender);
        confirmTransaction(transactionId);
        
        ownerActivityCount[msg.sender]++;
        lastActionTimestamp = block.timestamp;
    }
    
    function confirmTransaction(uint _transactionId)
        public
        onlyOwner
        transactionExists(_transactionId)
        notConfirmed(_transactionId)
        notExecuted(_transactionId)
        notCanceled(_transactionId)
        whenNotPaused
    {
        Transaction storage transaction = transactions[_transactionId];
        require(
            block.timestamp - transaction.submissionTime <= CONFIRMATION_TIMEOUT,
            "Confirmation timeout exceeded"
        );
        
        isConfirmed[_transactionId][msg.sender] = true;
        transaction.numConfirmations++;
        
        emit Confirmation(msg.sender, _transactionId, block.timestamp);
        
        ownerActivityCount[msg.sender]++;
        lastConfirmationTime[msg.sender] = block.timestamp;
        
        if (transaction.numConfirmations >= numConfirmationsRequired) {
            executeTransaction(_transactionId);
        }
    }
    
    function executeTransaction(uint _transactionId)
        public
        onlyOwner
        whenNotPaused
        notCanceled(_transactionId)
        transactionExists(_transactionId)
        notExecuted(_transactionId)
    {
        Transaction storage transaction = transactions[_transactionId];
        require(transaction.numConfirmations >= numConfirmationsRequired, "Not enough confirmations");
        
        if (transaction.value > 0) {
            uint256 today = block.timestamp / 1 days;
            if (today > lastDay) {
                lastDay = today;
                spentToday = 0;
            }
            
            require(spentToday + transaction.value <= dailyLimit, "Daily limit exceeded");
            spentToday += transaction.value;
        }
        
        transaction.executed = true;
        
        (bool success, bytes memory returnData) = transaction.to.call{value: transaction.value}(
            transaction.data
        );
        
        if (success) {
            emit Execution(_transactionId, block.timestamp);
            
            DailyOperation storage operation = dailyOperations[dailyOperationCount++];
            operation.value = transaction.value;
            operation.timestamp = block.timestamp;
            operation.initiator = msg.sender;
            operation.processed = true;
        } else {
            transaction.executed = false;
            emit ExecutionFailure(_transactionId, string(returnData));
        }
        
        lastActionTimestamp = block.timestamp;
    }
    
    function cancelTransaction(uint _transactionId)
        public
        onlyOwner
        transactionExists(_transactionId)
        notExecuted(_transactionId)
        notCanceled(_transactionId)
    {
        require(
            transactions[_transactionId].submitter == msg.sender,
            "Only submitter can cancel"
        );
        
        transactions[_transactionId].canceled = true;
        emit TransactionCanceled(_transactionId);
    }
    
    function addOwner(address _newOwner)
        public
        onlyOwner
        notNull(_newOwner)
        ownerDoesNotExist(_newOwner)
    {
        require(owners.length < MAX_OWNER_COUNT, "Max owners reached");
        isOwner[_newOwner] = true;
        owners.push(_newOwner);
        lastConfirmationTime[_newOwner] = block.timestamp;
        emit OwnerAddition(_newOwner);
    }
    
    function removeOwner(address _owner)
        public
        onlyOwner
        ownerExists(_owner)
    {
        require(owners.length > numConfirmationsRequired, "Cannot remove owner");
        isOwner[_owner] = false;
        for (uint i = 0; i < owners.length; i++) {
            if (owners[i] == _owner) {
                owners[i] = owners[owners.length - 1];
                owners.pop();
                break;
            }
        }
        emit OwnerRemoval(_owner);
    }
    
    function changeDailyLimit(uint256 _dailyLimit)
        public
        onlyOwner
    {
        dailyLimit = _dailyLimit;
        emit DailyLimitChange(_dailyLimit);
    }
    
    function togglePause()
        public
        onlyOwner
    {
        paused = !paused;
        if (paused) {
            emit ContractPaused(msg.sender);
        } else {
            emit ContractUnpaused(msg.sender);
        }
    }
    
    function getOwners()
        public
        view
        returns (address[] memory)
    {
        return owners;
    }
    
    function getTransactionCount()
        public
        view
        returns (uint)
    {
        return transactionCount;
    }
    
    function getTransaction(uint _transactionId)
        public
        view
        returns (
            address to,
            uint value,
            bytes memory data,
            bool executed,
            uint numConfirmations,
            uint256 submissionTime,
            address submitter,
            bool canceled,
            string memory description
        )
    {
        Transaction storage transaction = transactions[_transactionId];
        return (
            transaction.to,
            transaction.value,
            transaction.data,
            transaction.executed,
            transaction.numConfirmations,
            transaction.submissionTime,
            transaction.submitter,
            transaction.canceled,
            transaction.description
        );
    }
    
    function getConfirmationStatus(uint _transactionId, address _owner)
        public
        view
        returns (bool)
    {
        return isConfirmed[_transactionId][_owner];
    }
    
    function getBalance()
        public
        view
        returns (uint256)
    {
        return address(this).balance;
    }
    
    function getDailyStats()
        public
        view
        returns (
            uint256 _dailyLimit,
            uint256 _spentToday,
            uint256 _remainingToday
        )
    {
        uint256 today = block.timestamp / 1 days;
        if (today > lastDay) {
            return (dailyLimit, 0, dailyLimit);
        }
        return (dailyLimit, spentToday, dailyLimit - spentToday);
    }
}